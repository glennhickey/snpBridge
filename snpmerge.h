/*
 * Copyright (C) 2015 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.cactus
 */

#ifndef _SNPMERGE_H
#define _SNPMERGE_H

#include <string>
#include <vector>
#include <limits>
#include <map>
#include <stdexcept>

#include <sstream>

#include "vg/vg.hpp"
#include "Variant.h"

/** 
Let's say we have two adjacent snps, along with phasing information. 
The following haplotypes (for the pair of snps) are possible:
ref-ref
ref-alt
alt-ref
alt-alt
Chances are, the vcf file will only contain a subset of these posibilities.
But the vg graph constructed from the vcf will contain them all.  
So here we take as input a vg graph, and the vcf it was constructed from.  
Using the vcf we look at pairs of snps, and remove paths that don't 
correspond to haplotypes.  This is done with two operations
(credit: Benedict Paten):

1) 

Input: (bottom snps are reference)
A       G
 ACAFATA
T       C

Only ref-ref and alt-alt:
We duplicate reference between the snps,
and use it to make two parallel paths:

A-ACAFATA-G
 
T-ACAFATA-C
********************************88

2)

Only ref-ref annd alt-ref and ref-alt:

A-ACAFATA G
         X
T-ACAFATA-C

Duplicate reference as above.  This time
the dublicated reference segment can only 
connect to the reference base (C, lower-right),
but the original reference segment can connect
to both alts (C, G)

*/

class SNPMerge
{
public:

   SNPMerge();
   ~SNPMerge();
};

#endif
