/*
 * Copyright (C) 2015 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.cactus
 */

#include "snpbridge.h"

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
    // skip ahead until var2 doesn't overlap var1 or anything between
    bool breakOut = false;
    int prev_position = var1.position + var1.alleles[0].size();
    while (!breakOut && var2.position < prev_position)
    {
      cerr << "Skipping variant at " << var2.position << " because it "
           << "overlaps previous variant at position " << var1.position << endl;
      prev_position = max(prev_position,
                          (int)(var2.position + var2.alleles[0].size()));
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
#ifdef DEBUG
    cerr << "Linkcounts: ";
    for (int i = 0; i < _linkCounts.size(); ++i)
    {
      for (int j = 0; j < _linkCounts[i].size(); ++j)
      {
        cerr << "(" << i <<"-" << j << "=" << _linkCounts[i][j] << ") ";
      }
    }
    cerr << endl;
#endif

    for (int a1 = 1; a1 < var1.alleles.size(); ++a1)
    {
      for (int a2 = 1; a2 < var2.alleles.size(); ++a2)
      {
        // note can probably get what we need by calling once instead
        // of in loop....
        Phase phase = phaseRelation(var1, a1, var2, a2);

        if (phase != GT_OTHER)
        {
          makeBridge(a1, a2, phase);
          // we can get away with breaking here (and below) because results
          // mutually exclusive (see simplifying assumption in
          // phaseRelation()).  So as soon as we see a GT_AND or
          // GT_XOR, then everything else must be GT_OTHER
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

void SNPBridge::makeBridge(int allele1, int allele2, Phase phase)
{
#ifdef DEBUG
  cerr << allele1 << " " << phase2str(phase) << " " << allele2 << " detected at "
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
  Node* ref1 = _gv1.getGraphAllele(0).back();
  Node* node2 = _gv2.getGraphAllele(allele2).front();
  Node* ref2 = _gv2.getGraphAllele(0).front();

  // note we don't use references here because they get altered by
  // calls to create and destroy.
  vector<pair<int64_t, bool> > outEdges1 = _vg->edges_on_end[node1->id()];
  vector<pair<int64_t, bool> > inEdges2 = _vg->edges_on_start[node2->id()];

  // make sure there's no other way out of node1 but the
  // new bridges that we'll add
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

  // make sure there's no other way into node 2 than
  // the bridges we add
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

  // find the path between the two variant alleles along the
  // reference.  since we only deal with consecutive variants,
  // it's sufficient to stick this path between
  list<Node*> refPath;
  _gv1.getReferencePathTo(_gv2, refPath);

  // if there's no path, we assume the variants are directly adjacent
  // and just stick edges between them
  if (refPath.empty())
  {
    if (phase == GT_AND || phase == GT_FROM_REF || phase == GT_TO_REF)
    {
      _vg->create_edge(node1, node2, false, false);
#ifdef DEBUG
      cerr << "create adj alt-alt " << node1->id() << ", " << node2->id() << endl;
#endif
    }
    if (phase == GT_FROM_REF || phase == GT_XOR)
    {
      _vg->create_edge(ref1, node2, false, false);
#ifdef DEBUG
      cerr << "create adj ref-alt " << ref1->id() << ", " << node2->id() << endl;
#endif
    }
    if (phase == GT_TO_REF || phase == GT_XOR)
    {
      _vg->create_edge(node1, ref2, false, false);
#ifdef DEBUG
      cerr << "create adj alt-ref " << node1->id() << ", " << ref2->id() << endl;
#endif
    }
  }
  
  // otherwise, make a copy of ref path and stick that in between
  else
  {
    Node* prev = node1;
    Node* refPrev = ref1;
    for (auto refNode : refPath)
    {
      Node* cpyNode = _vg->create_node(refNode->sequence());
      _vg->create_edge(prev, cpyNode, false, false);
#ifdef DEBUG
      cerr << "create " << cpyNode->id() << endl;
      cerr << "create " << prev->id() << " -> " << cpyNode->id() << endl;
#endif
      prev = cpyNode;
      refPrev = refNode;
    }

    if (phase == GT_AND || phase == GT_FROM_REF || phase == GT_TO_REF)
    {
      _vg->create_edge(prev, node2, false, false);
#ifdef DEBUG
      cerr << "create alt-alt " << prev->id() << ", " << node2->id() << endl;
#endif
    }
    if (phase == GT_FROM_REF || phase == GT_XOR)
    {
      _vg->create_edge(refPrev, node2, false, false);
#ifdef DEBUG
      cerr << "create ref-alt " << refPrev->id() << ", " << node2->id() << endl;
#endif
    }
    if (phase == GT_TO_REF || phase == GT_XOR)
    {
      _vg->create_edge(prev, ref2, false, false);
#ifdef DEBUG
      cerr << "create alt-ref " << prev->id() << ", " << ref2->id() << endl;
#endif
    }
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
  
  bool to_ref = _linkCounts[allele1][0] > 0;
  bool from_ref = _linkCounts[0][allele2] > 0;
  bool to_alt = _linkCounts[allele1][allele2] > 0;

  bool to_other_alt = false;
  for (int i = 1; i < v2.alleles.size(); ++i)
  {
    if (i != allele2 && _linkCounts[allele1][i] > 0)
    {
      to_other_alt = true;
      break;
    }
  }
  
  bool from_other_alt = false;
  for (int i = 1; i < v1.alleles.size(); ++i)
  {
    if (i != allele1 && _linkCounts[i][allele2] > 0)
    {
      from_other_alt = true;
      break;
    }
  }

  // don't handle multi allele cases
  if (from_other_alt || to_other_alt)
  {
    return GT_OTHER;
  }
  
  if (to_alt)
  {
    if (!from_ref && !to_ref)
    {
      return GT_AND;
    }
    if (from_ref && !to_ref)
    {
      return GT_FROM_REF;
    }
    if (!from_ref && to_ref)
    {
      return GT_TO_REF;
    }
  }
  else
  {
    if (!to_ref)
    {
      cerr << "allele 1 " << allele1 << " to ref " << _linkCounts[allele1][0] << " and "
           << "to_ref " << to_ref << endl;
      cerr << "Alternate allele " << allele1 << " never seen in GT for "
           << "variant " << v1 << endl;
    }
    if (!from_ref)
    {
      cerr << "allele 2 " << allele2 << " from ref " << _linkCounts[0][allele2] << " and "
           << "from_ref " << from_ref << endl;
      cerr << "Alternate allele " << allele2 << " never seen in GT for "
           << "variant " << v2 << endl;
    }
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
    if (mapping.position().is_reverse() == true)
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
