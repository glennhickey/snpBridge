/*
 * Copyright (C) 2015 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.cactus
 */

#ifndef _SNPBRIDGE_H
#define _SNPBRIDGE_H

#include <string>
#include <vector>
#include <limits>
#include <map>
#include <stdexcept>
#include <sstream>

#include "vg/vg.hpp"
#include "Variant.h"
#include "graphvariant.h"

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

class SNPBridge
{
public:

   enum Phase {GT_AND, GT_XOR, GT_OTHER};

   SNPBridge();
   ~SNPBridge();

   /** iterate through adjacent snps and do merging, in place, in the 
    * vg graph */
   void processGraph(vg::VG* vg, vcflib::VariantCallFile* vcf, int offset,
                     int windowSize);


protected:

   /** Make a direct bridge from end of allele1 (in graph) to 
    * start of allele2.  If allele2 is not reference (0), 
    * make sure the new bridge is the only point of entry */
   void makeBridge(int allele1, int allele2);
   
   /** get the phasing relationship using the GT fields between
    * twp variants for a pair of two alternate alleles (>0).
    * GT_AND: both alleles present in all (non-ref) haplotypes
    * GT_XOR: exactly one allele present in all (non-ref) haplotypes
    * OTHER: other (no phasing information we can use to uncollapse).
    *
    * Important: need to call computeLinkCounts first. 
    */
   Phase phaseRelation(vcflib::Variant& v1, int allele1,
                       vcflib::Variant& v2, int allel2) const;


   /** Count the number of samples that have each pair of allele variants
    * on same haplotype, storing results in _linkCoutns array */
   void computeLinkCounts(vcflib::Variant& v1,
                          vcflib::Variant& v2);

   /** Resize and set matrix to 0 */
   void initLinkCounts(vcflib::Variant& v1,
                       vcflib::Variant& v2);

   /** Check length of reference path */
   int vgRefLength(vcflib::Variant& var) const;
   
protected:

   vg::VG* _vg;
   GraphVariant _gv1;
   GraphVariant _gv2;

   /** store the number of samples that have a pair variants on
    * the same allele.  These numbers can be used to tell if
    * to snps, for instance, area lways ref-ref / alt-alt.
    * so _linkCoutnos[1][1] = X would mean X samples
    * have alt1-alt1 for the two variants in consideration
    * on same chromosome. */
   std::vector<std::vector<int> > _linkCounts;
};

#endif
