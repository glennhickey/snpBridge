# snpBridge
Bridge adjacent snps in a variant graph constructed from vcf, using genotype information in said vcf.  The goal is to remove paths from the graph that are not present in any haplotypes implied by the VCF.

## Limitations
* **Will only work on vg files created with vg construct -f -v**
* Only adjacent pairs of VCF variants are considered
* Overlapping variants in VCF will be skipped
* Coordinate of first vg position must be known (and passed with -o).  This would normall be the start position given in vg construct -R to create the initial graph...
* Multiallelic variants are handled, but not combinatorially.  Ie only the two cases above are considered for each allele.  

## Method

Adjacent SNPs are bridged if they are within the window size and one of the following two cases apply:

### Case 1: Linked Alternate Alleles

If given the genotype information, every sample allele that has alternate at SNP 1 also has alternate at SNP 2, then the graph is changed so that any path going through the alternate allele at SNP1 must also go through the alternate allele at SNP 2.  

#### Example: Before
![altalt_orig](https://raw.githubusercontent.com/glennhickey/snpBridge/development/doc/altalt_orig.png)
#### Example: After
![altalt_orig](https://raw.githubusercontent.com/glennhickey/snpBridge/development/doc/altalt_bridge.png)

### Case 2: Linked Alternate and Reference Alleles

If given the genotype information, every sample allele that has alternate at SNP 1 has reference at SNP 2, then the graph is changed so that any path going through the alternate allele at SNP1 must also go through the reference allele at SNP 2.  

#### Example: Before
![altalt_orig](https://raw.githubusercontent.com/glennhickey/snpBridge/development/doc/altref_orig.png)
#### Example: After
![altalt_orig](https://raw.githubusercontent.com/glennhickey/snpBridge/development/doc/altref_bridge.png)


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
