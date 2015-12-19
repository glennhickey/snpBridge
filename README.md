# snpBridge
Bridge adjacent snps in a variant graph constructed from vcf, using genotype information in said vcf.  The goal is to remove paths from the graph that are not present in any haplotypes implied by the VCF.

1. If two alternate alleles in consecutive variants are linked in all samples, they will be bridged so that all paths in the graph will either contain both or none.
2. If an alternate allele is linked to the *reference* at the next variant in all samples, it will be bridged so that all paths in the graph that contain it never contain any alts at the next variant.


## Limitations
* Will only work on vg files created with vg construct -f -v
* Only adjacent pairs of VCF variants are considered
* Overlapping variants in VCF will be skipped
* Coordinate of first vg position must be known (and passed with -o).  This would normall be the start position given in vg construct -R to create the initial graph...
* Multiallelic variants are handled, but not combinatorially.  Ie only the two cases above are considered for each allele.  
