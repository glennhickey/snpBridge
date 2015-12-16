/*
 * Copyright (C) 2015 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.cactus
 */

#ifndef _GRAPHVARIANT_H
#define _GRAPHVARIANT_H

#include <string>
#include <vector>
#include <limits>
#include <map>
#include <stdexcept>

#include <sstream>

#include "vg/vg.hpp"
#include "Variant.h"

/** 
Maintain a mapping between a vcf variant and the vg graph
positions to which each of its alleles correspond. 

Note: This only works on vg files created with vg construct -f
Otherwise multibase events will trigger errors.  Maybe add
option to skip these...
*/

class GraphVariant
{
public:

   enum Cat {SNP, DEL, INS, INDEL, REFONLY};
   
   GraphVariant();
   GraphVariant(vcflib::Variant& var);
   
   ~GraphVariant();

   /** reset path and set offset */
   void init(int offset);
   
   /** scan forward along the path until we hit given Variant.
    */
   void loadVariant(vg::VG* vg, vcflib::Variant& var);

   /** how many alleles, reference included, at current variant 
    */
   int getNumAlleles() const;

   /** get vcf allele by number */
   const string& getVCFAllele(int i) const;

   /** get vg allele by number */
   const std::list<vg::Node*>& getGraphAllele(int i) const;

   /** what kind of variant.  */
   Cat varCat(vcflib::Variant& var) const;

   /** compare substring [o1, o1+len) of s1 to s2 beginning at o2.
    * case insensitive.  if either of the string not long enough
    * returns false*/
   static bool istreq(const std::string& s1, const std::string& s2,
                      int o1 = 0, int o2 = 0, int len = -1);

protected:

   /** read in the vg alleles by searching graph bubbles */
   void loadAlleles();
   
protected:
   
   vg::VG* _vg;
   std::list<vg::Mapping>* _path;
   std::list<vg::Mapping>::const_iterator _mappingIt;
   int _refIdx;
   vcflib::VariantCallFile* _vcf;
   vcflib::Variant _var;
   int _offset;
   Cat _cat;

   std::vector<std::list<vg::Node*> > _graphAlleles;
   
};

#endif
