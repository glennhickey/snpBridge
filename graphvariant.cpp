/*
 * Copyright (C) 2015 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.cactus
 */
#include "graphvariant.h"

#include "vg/vg.hpp"
#include "Variant.h"

using namespace vcflib;
using namespace vg;
using namespace std;

GraphVariant::GraphVariant() : _vg(NULL), _path(NULL), _refIdx(-1),
                               _vcf(NULL), _offset(0), _cat(REFONLY)
{
}

GraphVariant::~GraphVariant()
{
}

void GraphVariant::init(int offset)
{
  _vg = NULL;
  _path = NULL;
  _refIdx = -1;
  _vcf = NULL;
  _offset = offset;
}

void GraphVariant::loadVariant(VG* vg, Variant& var)
{
  _vg = vg;
  _var = var;
  _cat = varCat(var);
  
  // no path, let's find it and start at the beginning
  if (_path == NULL || _refIdx < 0)
  {
    if (_vg->paths.has_path(var.sequenceName) == false)
    {
      stringstream ss;
      ss << "Unable to find path for " << var.sequenceName << " in vg file";
      throw runtime_error(ss.str());
    }
    _path = &_vg->paths.get_path(var.sequenceName);
    _mappingIt = _path->begin();
    _refIdx = 0;
  }

  // scan forward in vg along path until are variant reference falls within
  // a mapping
  int mappingLen = -1;
  int vcfPos = var.position - _offset;
  for (; _mappingIt != _path->end(); ++_mappingIt, _refIdx += mappingLen)
  {
    if (_mappingIt->is_reverse() == true)
    {
      throw(runtime_error("Reverse Mapping not supported"));
    }
    if (_mappingIt->edit_size() > 1 || (
          _mappingIt->edit_size() == 1 && _mappingIt->edit(0).from_length() !=
          _mappingIt->edit(0).to_length()))
    {
      stringstream ss;
      ss << pb2json(*_mappingIt) << ": Only mappings with a single trvial edit"
         << " supported in ref path";
      throw runtime_error(ss.str());
    }

    mappingLen =
       _vg->get_node(_mappingIt->position().node_id())->sequence().length();
  
    if (_refIdx + mappingLen > vcfPos)
    {
      break;
    }
  }

  if (mappingLen < 1 || vcfPos < _refIdx || vcfPos >= _refIdx + mappingLen)
  {
    stringstream ss;
    ss << "Variant at position " << var.sequenceName << ":" << var.position
       << " not found in graph";
    throw runtime_error(ss.str());
  }

  // ok, now we should have found the exact place in the vg.  let's load
  // the alleles
  loadAlleles();
}

void GraphVariant::loadAlleles()
{
#ifdef DEBUG
  cerr << "1st vg ref node " << _mappingIt->position().node_id()
       << " _var position " << _var.position << " ref "
       << _var.alleles[0] << endl;
#endif
  // this is the structure we'll fill.  ith item corresponds to
  // path in graph that matches ith allele
  _graphAlleles.resize(1);
  _graphAlleles[0].clear();
  
  // find path in the graph corresponding to reference allele
  // because we assume vg construct -f used, the allele shouold
  // be exactly represented by a path (with no offsets)
  string vgRefPath;
  std::list<vg::Mapping>::const_iterator it = _mappingIt;
  for (; vgRefPath.length() < _var.alleles[0].length() &&
          it != _path->end(); ++it)
  {
    Node* node =  _vg->get_node(it->position().node_id());
    vgRefPath += node->sequence();
    _graphAlleles[0].push_back(node);
  }
  if (!istreq(vgRefPath, _var.alleles[0]))
  {
    stringstream ss;
    ss << "vg path beginning at node " << _graphAlleles[0].front()->id()
       << " has sequence " << vgRefPath << " which does not "
       << " match vcf reference for this record: " << _var;
    throw runtime_error(ss.str());
  }

#ifdef DEBUG
  cerr << "VCF: " << _var.sequenceName << "\t" << (_var.position - _offset)
       << "\t" << _var.alleles[0] << " VG:";
  for (auto& n : _graphAlleles[0])
  {
    cerr << "\t" << pb2json(*n);
  }
  cerr << endl;
#endif
  
  // now, our variants will be in the set of siblings
  // (nodes that share neighbouring sides on both ends
  // as our reference)  
  set<NodeTraversal> toSibs = _vg->full_siblings_to(_graphAlleles[0].front());
  set<NodeTraversal> fromSibs = _vg->full_siblings_from(_graphAlleles[0].back());
  set<NodeTraversal> sibs;

  //set_intersection(toSibs.begin(), toSibs.end(),
  //                 fromSibs.begin(), fromSibs.end(), sibs.begin());
  // Weird: above won't compile.  just brute force it
  for (auto& nt : toSibs)
  {
    if (fromSibs.find(nt) != fromSibs.end())
    {
      sibs.insert(nt);
    }
  }
      
  // search siblings for remaining alleles.  expecting exact mathc
  // of vg node to vcf allele
  for (int i = 1; i < _var.alleles.size(); ++i)
  {
    for (auto& nt : sibs)
    {
      if (istreq(_var.alleles[i], nt.node->sequence()))
      {
        _graphAlleles.push_back(list<Node*>(1, nt.node));
      }
    }
    
    // not finding a node for this allele is an error
    if (_graphAlleles.size() != i + 1)
    {
      stringstream ss;
      ss << "Unable to find vg node for allele " << i << " of " <<_var;
      throw runtime_error(ss.str());
    }
  }
}

GraphVariant::Cat GraphVariant::varCat(Variant& var) const
{
  int ref_len = var.alleles[0].length();
  bool ins = false;
  bool del = false;
  bool snp = false;
  for (int i = 1; i < var.alleles.size(); ++i)
  {
    if (var.alleles[i].length() == ref_len)
    {
      snp = true;
    }
    else if (var.alleles[i].length() > ref_len)
    {
      ins = true;
    }
    else if (var.alleles[i].length() < ref_len)
    {
      del = true;
    }
  }

  if (ins)
  {
    return del ? INDEL : INS;
  }
  else if (del)
  {
    return ins ? INDEL : DEL;
  }
  else if (snp)
  {
    return SNP;
  }
  
  return REFONLY;
}

int GraphVariant::getNumAlleles() const
{
  return _var.alleles.size();
}

const string& GraphVariant::getVCFAllele(int i) const
{
  return _var.alleles[i];
}

const list<Node*>& GraphVariant::getGraphAllele(int i) const
{
  return _graphAlleles[i];
}

const Variant& GraphVariant::getVariant() const
{
  return _var;
}

void GraphVariant::getReferencePathTo(const GraphVariant& other,
                                      list<Node*>& outPath) const
{
  outPath.clear();

  // scan to end of reference allele
  std::list<vg::Mapping>::const_iterator i = _mappingIt;
  for (; i != _path->end() &&
          i->position().node_id() != _graphAlleles[0].back()->id(); ++i);
  assert(i != _path->end());
  assert(i->position().node_id() == _graphAlleles[0].back()->id());

  // our destination:
  Node* otherHead = other._graphAlleles[0].front();

  // walk one step, then start keeping track of our path
  // (all nodes between _graphAlleles[0].back() and
  //  other._graphAlleles[0].front(), exclusive)
  for (++i; i != _path->end() && i->position().node_id() != otherHead->id(); ++i)
  {
    outPath.push_back(_vg->get_node(i->position().node_id()));
  }

  // no way we should get to end since other is between it and us
  assert(i != _path->end());
}

bool GraphVariant::overlaps(const GraphVariant& other) const
{
  return _var.position + _var.alleles[0].size() > other._var.position;
}

bool GraphVariant::istreq(const string& s1, const string& s2,
                          int o1, int o2, int len)
{
  if (len < 0)
  {
    len = max(s1.length() - o1, s2.length() - o2);
  }
  assert(len > 0);
  for (int i = 0; i < len; ++i)
  {
    if (o1 + i >= s1.length() || o2 + i >= s2.length() ||
        toupper(s1[o1 + i]) != toupper(s2[o2 + i]))
    {
      return false;
    }
  }
  return true;
}

ostream& operator<<(ostream& os, const GraphVariant& gv)
{
  const Variant v = gv.getVariant();
  os << "GV:[" << v.sequenceName << ":" << v.position;
  for (int i = 0; i < v.alleles.size(); ++i)
  {
    os << "(" << i << ") vcf=" << v.alleles[i] << " | vg=";
    for (auto n : gv.getGraphAllele(i))
    {
      os << n->id() << ",";
    }
    os << "; ";
  }
  os << "]";
  return os;
}
