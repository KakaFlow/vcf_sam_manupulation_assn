## VCF file manupulation
#1.Structure of vcf file
bcftools view files/sample.vcf.gz

#2. Header content
bcftools view -h files/sample.vcf.gz

#3. How many samples are in the file.
bcftools query -l files/sample.vcf.gz | wc -l

#4. How many variants are in the file
bcftools query -f '%ALT\n' files/sample.vcf.gz | wc -l`

#5. How would you extract the chromosome, position, QualByDepth and RMSMappingQuality fields? Save the output to a tab-delimited file
bcftools query -f '%CHROM\t%POS\t%QUAL\t%MQ\n' files/sample.vcf.gz > out_dir/answer5.tsv

#6. Extract data that belongs to chromosomes 2,4 and MT
zcat files/sample.vcf.gz | awk -F "\t" '$1~/^(2|4|MT)$/{print $0}'

#7. Print out variants that do not belong to chr20:1-30000000
zcat files/sample.vcf.gz |  grep '^[^#]' | awk '!($1 == 20 && $2 >= 1 && $2 <= 30000000)'

#8. Extract variants that belong to SRR13107019
bcftools query -f '%ALT\n' -s SRR13107019 files/sample.vcf.gz

#9. Filter out variants with a QualByDepth above 7
bcftools query -f '[%ALT\t%QD\n]' files/sample.vcf.gz | awk '$2>7'

10. How many contigs are referred to in the file. Check the header section
11. Comment on the eighth and ninth columns of the file
12. Extract data on the read depth of called variants for sample SRR13107018
13. Extract data on the allele frequency of alternate alleles. Combine this data with the
chromosome and position of the alternate allele
