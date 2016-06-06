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

#include "vg/src/vg.hpp"
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
    correspond to haplotypes.  This is done by inserting a bridge with 
    appropriate connectors
    (credit: Benedict Paten):

    ex: 

    Input: (bottom snps are reference)
    A       G
    ACAFATA
    T       C

    We duplicate reference between the snps,
    and use it to make two parallel paths:

    A-ACAFATA-G
 
    T ACAFATA C

    Edges then filled in to represent linking in VCF...
    ex: If a haplotype contains the C in alt to iff if has the T at alt 1:
    (GT_AND)

    A-ACAFATA-G
 
    T-ACAFATA-C

    but if no haplotype has both the T and C alt:
    (GT_XOR)
    A-ACAFATA-G
             X
    T-ACAFATA C

    or no hapltoype has T-Alt and G-Ref
    (GT_FROM_REF)
    A-ACAFATA-G
             \
    T-ACAFATA-C

    or no hapltoype has A-Ref and T-Alt
    (GT_TO_REF)
    A-ACAFATA-G
             /
    T-ACAFATA-C
*/

class SNPBridge
{
public:

   enum Phase {GT_AND = 0, // only link from Alt1 to Alt2
               GT_XOR, // only link from Alt1 to ref and ref to Alt2
               GT_FROM_REF, // link from Alt1 to Alt2 and ref to Alt2
               GT_TO_REF, // link from Alt1 to Alt2 and Alt1 to ref
               GT_OTHER}; // all links (leave alone)

   SNPBridge();
   ~SNPBridge();

   /** iterate through adjacent snps and do merging, in place, in the 
    * vg graph */
   void processGraph(vg::VG* vg, vcflib::VariantCallFile* vcf, int offset,
                     int windowSize);


protected:

   /** Make a direct bridge from end of allele1 (in graph) to 
    * start of allele2. */
   void makeBridge(int allele1, int allele2, Phase phase);
   
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

inline std::string phase2str(SNPBridge::Phase phase)
{
  static string s[] = {"GT_AND", "GT_XOR", "GT_FROM_REF", "GT_TO_REF", "GT_OTHER"}; 
  return s[phase];
}
#endif
