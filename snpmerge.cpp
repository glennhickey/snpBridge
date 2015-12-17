/*
 * Copyright (C) 2015 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.cactus
 */

#include "snpmerge.h"

#include "vg/vg.hpp"
#include "Variant.h"

using namespace vcflib;
using namespace vg;
using namespace std;

SNPMerge::SNPMerge()
{
}

SNPMerge::~SNPMerge()
{
}

void SNPMerge::processGraph(VG* vg, VariantCallFile* vcf, int offset)
{
  _gv1.init(offset);
  _gv2.init(offset);
  
  Variant var1(*vcf);
  Variant var2(*vcf);
  if (!vcf->getNextVariant(var1))
  {
    // empty file
    cerr << "No variants found in VCF" << endl;
    return;
  }
  _gv1.loadVariant(vg, var1);

  while (vcf->getNextVariant(var2))
  {
    _gv2.loadVariant(vg, var2);

    cerr << "v1 " << var1.sequenceName << " " << var1.position
         << " v2 " << var2.sequenceName << " " << var2.position << endl;

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
          uncollapse_and(a1, a2);
        }
        else if (phase == GT_XOR)
        {
          uncollapse_xor(a1);
        }
        // we can get away with breaking here because results
        // mutually exclusive (see simplifying assumption in
        // phaseRelation()).  So as soon as we see a GT_AND or
        // GT_XOR, then everything else must be GT_OTHER
        else
        {
          break;
        }
      }
    }
    
    swap(var1, var2);
    swap(_gv1, _gv2);
  }  
}

void SNPMerge::uncollapse_and(int allele1, int allele2)
{
#ifdef DEBUG
  cerr << allele1 << " AND " << allele2 << " detected at "
       << _gv1.getVariant().position << endl;
#endif
  
}

void SNPMerge::uncollapse_xor(int allele1)
{
#ifdef DEBUG
  cerr << allele1 << " OR x" << " detected at "
       << _gv1.getVariant().position << endl;
#endif
  
}

SNPMerge::Phase SNPMerge::phaseRelation(Variant& v1, int allele1,
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
    cerr << "No phasing information found for var1=" << v1
         << endl << " and v2 " << v2 << endl;
    return GT_OTHER;
  }
  
  // so exclusively to allele 2 will be a GT_AND
  if (toAllele2 == total)
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

void SNPMerge::computeLinkCounts(Variant& v1, Variant& v2)
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

void SNPMerge::initLinkCounts(Variant& v1,
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
