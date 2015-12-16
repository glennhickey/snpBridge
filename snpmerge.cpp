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
    
    swap(var1, var2);
    swap(_gv1, _gv2);
  }  
}
