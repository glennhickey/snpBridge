/*
 * Copyright (C) 2015 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.cactus
 */

#include "snpbridge.h"

#include "vg/vg.hpp"
#include "Variant.h"

using namespace vcflib;
using namespace vg;
using namespace std;

SNPBridge::SNPBridge() : _vg(NULL)
{
}

SNPBridge::~SNPBridge()
{
}

void SNPBridge::processGraph(VG* vg, VariantCallFile* vcf, int offset,
                            int windowSize)
{
  _vg = vg;
  _gv1.init(offset);
  _gv2.init(offset);
  
  Variant var1(*vcf);
  Variant var2(*vcf);
  // skip to first variant after offset
  for (int vcfPos = -1; vcfPos < offset; vcfPos = var1.position)
  {
    if (!vcf->getNextVariant(var1))
    {
      // empty file
      cerr << "No variants found in VCF" << endl;
      return;
    }
  }
  _gv1.loadVariant(vg, var1);

  int graphLen = vgRefLength(var1);


  for (; vcf->getNextVariant(var2); swap(var1, var2), swap(_gv1, _gv2))
  {
    // skip ahead until we're not overlapping
    bool breakOut = false;
    while (!breakOut && var2.position < var1.position + var1.alleles[0].size())
    {
      cerr << "Skipping variant at " << var2.position << " because it "
           << "overlaps previous variant at position " << var1.position << endl;
      breakOut = !vcf->getNextVariant(var2);
    }
    if (breakOut)
    {
      break;
    }

    if (var2.position >= offset + graphLen)
    {
      // stop after end of vg
      break;
    }
    
    _gv2.loadVariant(vg, var2);
    
    if (var2.position - (var1.position + var1.alleles[0].length() - 1) >
        windowSize)
    {
      // skip because further than window size
      continue;
    }

#ifdef DEBUG
    cerr << "\nv1 " << _gv1 << endl << "v2 " << _gv2 << endl;
#endif

    computeLinkCounts(var1, var2);

    for (int a1 = 1; a1 < var1.alleles.size(); ++a1)
    {
      for (int a2 = 1; a2 < var2.alleles.size(); ++a2)
      {
        // note can probably get what we need by calling once instead
        // of in loop....
        Phase phase = phaseRelation(var1, a1, var2, a2);

        if (phase == GT_AND)
        {
          makeBridge(a1, a2);
          // we can get away with breaking here (and below) because results
          // mutually exclusive (see simplifying assumption in
          // phaseRelation()).  So as soon as we see a GT_AND or
          // GT_XOR, then everything else must be GT_OTHER
          break;
        }
        else if (phase == GT_XOR)
        {
          makeBridge(a1, 0);
          break;
        }
        else
        {
#ifdef DEBUG
          cerr << a1 << " OTHER " << a2 << " detected at "
               << _gv1.getVariant().position << " ";
          for (int i = 0; i < _linkCounts.size(); ++i)
          {
            cerr << "(";
            for (int j = 0; j < _linkCounts[i].size(); ++j)
            {
              cerr << _linkCounts[i][j] << ",";
            }
            cerr << ") ";
          }
          cerr << endl;
#endif
        }
      }
    }
  }
}

void SNPBridge::makeBridge(int allele1, int allele2)
{
#ifdef DEBUG
  cerr << allele1 << " UNC " << allele2 << " detected at "
       << _gv1.getVariant().position << " ";
  for (int i = 0; i < _linkCounts.size(); ++i)
  {
     cerr << "(";
     for (int j = 0; j < _linkCounts[i].size(); ++j)
     {
       cerr << _linkCounts[i][j] << ",";
     }
     cerr << ") ";
  }
  cerr << endl;
#endif

  Node* node1 = _gv1.getGraphAllele(allele1).back();
  Node* node2 = _gv2.getGraphAllele(allele2).front();

  // note we don't use references here because they get altered by
  // calls to create and destroy.
  vector<pair<int64_t, bool> > outEdges1 = _vg->edges_on_end[node1->id()];
  vector<pair<int64_t, bool> > inEdges2 = _vg->edges_on_start[node2->id()];

  // make sure there's no other way out of node1 but the
  // new bridge that we'll add
  for (auto p : outEdges1)
  {
#ifdef DEBUG
    cerr << "destroy1 " << node1->id() << " " << true << ", " << p.first << " " << p.second << endl;
#endif
    Edge* edge = _vg->get_edge(NodeSide(node1->id(), true),
                               NodeSide(p.first, p.second));
    assert(edge != NULL);
    _vg->destroy_edge(edge);
  }
  
  // if dealting with an alt-alt case, make sure no other
  // way into second alt than the new bridge we'll add
  // (for alt-ref case, just delete any edge from prev alt)
  if (allele2 != 0)
  {
    for (auto p : inEdges2)
    {
#ifdef DEBUG
      cerr << "destroy2 " << p.first << " " << !p.second << ", " << node2->id() << " " << false << endl;
#endif
      Edge* edge = _vg->get_edge(NodeSide(p.first, !p.second),
                                 NodeSide(node2->id(), false));
      if (p.first == node1->id())
      {
        // should have been deleted above
        assert(edge == NULL);
      }
      else
      {
        assert(edge != NULL);
        _vg->destroy_edge(edge);
      }
    }
  }

  // find the path between the two variant alleles along the
  // reference.  since we only deal with consecutive variants,
  // it's sufficient to stick this path between
  list<Node*> refPath;
  _gv1.getReferencePathTo(_gv2, refPath);

  // if there's no path, we assume the variants are directly adjacent
  // and just stick and edge between them
  if (refPath.empty())
  {
    _vg->create_edge(node1, node2, false, false);
#ifdef DEBUG
    cerr << "create " << node1->id() << ", " << node2->id() << endl;
#endif
  }
  // otherwise, make a copy of ref path and stick that in between
  else
  {
    Node* prev = node1;
    for (auto refNode : refPath)
    {
      Node* cpyNode = _vg->create_node(refNode->sequence());
      _vg->create_edge(prev, cpyNode, false, false);
#ifdef DEBUG
      cerr << "create " << cpyNode->id() << endl;
      cerr << "create " << prev->id() << " -> " << cpyNode->id() << endl;
#endif
      prev = cpyNode;
    }
    _vg->create_edge(prev, node2, false, false);
#ifdef DEBUG
    cerr << "create " << prev->id() << ", " << node2->id() << endl;
#endif
  }
}

SNPBridge::Phase SNPBridge::phaseRelation(Variant& v1, int allele1,
                                        Variant& v2, int allele2) const
{
  // this is where we could take into account allele
  // frequencies to, for example, ignore really rare alleles.
  // But for now, we only do an all or nothing -- ie
  // the variants are alt-alt only if there isn't a single sample
  // saying otherwise.

  // total number of haplos with a1 that also have a2
  int toAllele2 = _linkCounts[allele1][allele2];;
  // total number of haplos with a1 that go to ref
  int toRef = _linkCounts[allele1][0];
  // total number of implied haplotypes that contain allele1
  int total = 0;
  for (int i = 0; i < v2.alleles.size(); ++i)
  {
    total += _linkCounts[allele1][i];
  }

  if (total < 1)
  {
    cerr << "Alternate allele " << allele1 << " never seen in GT for "
         << "variant " << v1 << endl;
    // the best we can do is cut path to next alternate.
    // though really we want to remove such variants
    // entirely
    return GT_XOR;
  }
  
  // so exclusively to allele 2 will be a GT_AND
  else if (toAllele2 == total)
  {
    return GT_AND;
  }
  // and exclusive to ref with be a GT_XOR
  // (Note a simplifiying assumption here in the
  //  presence of multiple alts.  Don't consider
  //  cases like never goes to alt2 but can go
  //  to alt1 or reference)
  else if (toRef == total)
  {
    return GT_XOR;
  }
  
  return GT_OTHER;
}

void SNPBridge::computeLinkCounts(Variant& v1, Variant& v2)
{
  // make our matrix and set to 0
  initLinkCounts(v1, v2);
  
  assert(!v1.alleles.empty() && !v2.alleles.empty());

  for (auto& sample : v1.sampleNames)
  {
    // treat missing GT information in one variant with respect
    // to the other as a warning. But count all possible links
    // once so it will never get phased. 
    if (v2.samples.find(sample) == v2.samples.end())
    {
      cerr << "Warning: Sample " << sample << " not found in variant " << v2
           << ". Assuming unphased" << endl;
      for (auto& i : _linkCounts)
      {
        for (auto& j : i)
        {
          ++j;
        }
      }
    
      continue;
    }  

    // parse out GT info for sample (it'll look like 0|1 etc.)
    string& gt1 = v1.samples[sample]["GT"].front();
    string& gt2 = v2.samples[sample]["GT"].front();
    vector<string> gtspec1 = split(gt1, "|");
    vector<string> gtspec2 = split(gt2, "|");

    // don't expect this but check in case
    if (gtspec1.size() != gtspec2.size())
    {
      stringstream ss;
      ss << "Sample " << sample << "has different ploidy in "
         << v1 << " and " << v2 << ". This is not supported";
      throw runtime_error(ss.str());
    }

    // so if the GT info looks like 0|1 1|0 2|1 etc.  we have
    // two "chromosomes", and iterate each.  we count a link
    // looking at the two variants, and recording what we see
    // at the same chromosome at the same sample. 
    for (int chrom = 0; chrom < gtspec1.size(); ++chrom)
    {
      // treat . as wildcard
      if (gtspec1[chrom] == "." && gtspec2[chrom] == ".")
      {
        // two .'s mean we see everything
        for (int g1 = 0; g1 < v1.alleles.size(); ++g1)
        {
          for (int g2 = 0; g2 < v2.alleles.size(); ++g2)
          {
            ++_linkCounts[g1][g2];
          }
        }
      }
      else if (gtspec1[chrom] == ".")
      {
        // g1 == . means we see all combindations of g1 with
        // given value of g2
        int g2;
        convert(gtspec2[chrom], g2);
        for (int g1 = 0; g1 < v1.alleles.size(); ++g1)
        {
          ++_linkCounts[g1][g2];
        }
      }
      else if (gtspec2[chrom] == ".")
      {
        // g2 == . means we see all combinations of g2 with
        // given value of g1
        int g1;
        convert(gtspec1[chrom], g1);
        for (int g2 = 0; g2 < v2.alleles.size(); ++g2)
        {
          ++_linkCounts[g1][g2];
        }
      }
      else
      {
        // normal case.  update the link count for the indexes
        // found on the give allele.  Example
        // Sample NA12878 has GT 0|1 for var1 and 0|0 for var2
        // then we update _linkCounts[0][1] + 1 (when chrom = 0)
        // and _linkCounots[0][0] + 1 (when chrom = 1)
        int g1, g2;
        convert(gtspec1[chrom], g1);
        convert(gtspec2[chrom], g2);
        ++_linkCounts[g1][g2];
      }
    }
  }
}

void SNPBridge::initLinkCounts(Variant& v1,
                              Variant& v2)
{
  // set _linkCounts to 0 and make sure it's at least
  // big enough to hold our alleles
  for (int i = 0; i < v1.alleles.size(); ++i)
  {
    if (_linkCounts.size() < v1.alleles.size())
    {
      _linkCounts.push_back(vector<int>(v2.alleles.size(), 0));
    }
    else
    {
      if (_linkCounts[i].size() < v2.alleles.size())
      {
        _linkCounts[i].resize(v2.alleles.size());
      }
      _linkCounts[i].assign(v2.alleles.size(), 0);
    }
  }  
}

int SNPBridge::vgRefLength(Variant& var) const
{
  // duplicating some code from the built in traversal of graphvariant,
  // but it's nice to have path length at outset to make scope checking
  // simpler (to relax assumption that vg contains whole vcf). 
  if (_vg->paths.has_path(var.sequenceName) == false)
  {
    stringstream ss;
    ss << "Unable to find path for " << var.sequenceName << " in vg file";
    throw runtime_error(ss.str());
  }
  int len = 0;
  list<Mapping>& path = _vg->paths.get_path(var.sequenceName);
  for (auto& mapping : path)
  {
    if (mapping.is_reverse() == true)
    {
      throw(runtime_error("Reverse Mapping not supported"));
    }
    if (mapping.edit_size() > 1 || (
          mapping.edit_size() == 1 && mapping.edit(0).from_length() !=
          mapping.edit(0).to_length()))
    {
      stringstream ss;
      ss << pb2json(mapping) << ": Only mappings with a single trvial edit"
         << " supported in ref path";
      throw runtime_error(ss.str());
    }
    len += _vg->get_node(mapping.position().node_id())->sequence().length();
  }
  return len;
}
