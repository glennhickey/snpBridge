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

## Usage

     snpBridge [options] VGFILE VCFFILE

**options**

    -h, --help          print this help message
    -w, --window-size N maximum distance between adjacent snps to be merged (default=50)
    -o, --offset N      vcf-coordinate of first position in vg path (default=1)

## Exmaple

These commands will process the first 500 bases of the BRCA1 region in GRCh38.  Need the relevant vcf and fasta file (chromosome 17).  The merged and original graphs will be drawn in PDF

     WIN=500
     START=43044345
     END=43044646

     vg construct -v brca1.vcf.gz -r 17.fa -R 17:${START}-${END} -f > test.vg
     gzip -d brca1.vcf.gz -c > test.vcf
     snpBridge test.vg test.vcf -o ${START} -w ${WIN} > merge.vg 2> log.txt
     snpBridge test.vg test.vcf -o ${START} -w 0 > orig.vg 2> log2.txt

     vg view -pd merge.vg | dot -Tpdf > merge.pdf
     vg view -pd orig.vg | dot -Tpdf > orig.pdf
