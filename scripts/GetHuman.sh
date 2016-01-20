#!/bin/sh
# =============================================================================
# GET GRC HUMAN GENOME
for((x=1;x!=23;++x)); do wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/Assembled_chromosomes/seq/hs_ref_GRCh38_chr$x.fa.gz ; done
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/Assembled_chromosomes/seq/hs_ref_GRCh38_chrX.fa.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/Assembled_chromosomes/seq/hs_ref_GRCh38_chrY.fa.gz
for((x=1;x!=23;++x)); do zcat hs_ref_GRCh38_chr$x.fa.gz | grep -v ">" | tr -d -c "ACGT" > GRC$x; done
zcat hs_ref_GRCh38_chrX.fa.gz | grep -v ">" | tr -d -c "ACGT" > HSC23 ;
zcat hs_ref_GRCh38_chrY.fa.gz | grep -v ">" | tr -d -c "ACGT" > HSC24 ;
cat HSC* > HS.acgt;
rm -f *.fa.gz ;
#rm -f GRC*
# =============================================================================
