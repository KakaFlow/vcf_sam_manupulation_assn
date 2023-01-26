## VCF file manipulation

#1.Structure of vcf file
zcat files/sample.vcf.gz
#Or
bcftools view files/sample.vcf.gz

#2. Header content
zcat files/sample.vcf.gz |  grep '^##'
#Or
bcftools view -h files/sample.vcf.gz

#3. How many samples are in the file.
zcat files/sample.vcf.gz | grep -m1 "^#CHROM" | cut -f 10- | tr "\t" "\n"  | wc -l
#Or
bcftools query -l files/sample.vcf.gz | wc -l

#4. How many variants are in the file
zcat files/sample.vcf.gz |  grep -c '^[^#]'
#Or
bcftools query -f '%ALT\n' files/sample.vcf.gz | wc -l`

#5. How would you extract the chromosome, position, QualByDepth and RMSMappingQuality fields? Save the output to a tab-delimited file
zcat files/sample.vcf.gz |  grep '^[^#]' | cut -f 1,2,6,8 | awk 'match($0, /MQ=([0-9]+\.[0-9]+)\;/){ $4=substr($0, RSTART, RLENGTH-1) }1' | awk '{ sub("MQ=","",$4); print }' OFS='\t' > out_dir/answer5.tsv`
#Or
bcftools query -f '%CHROM\t%POS\t%QUAL\t%MQ\n' files/sample.vcf.gz > out_dir/answer5.tsv

#6. Extract data that belongs to chromosomes 2,4 and MT
zcat files/sample.vcf.gz | awk -F "\t" '$1~/^(2|4|MT)$/{print $0}'

#7. Print out variants that do not belong to chr20:1-30000000
zcat files/sample.vcf.gz |  grep '^[^#]' | awk '!($1 == 20 && $2 >= 1 && $2 <= 30000000)'

#8. Extract variants that belong to SRR13107019
bcftools query -f '%ALT\n' -s SRR13107019 files/sample.vcf.gz

#9. Filter out variants with a QualByDepth above 7
zcat files/sample.vcf.gz |  grep '^[^#]' | awk 'match($0, /QD=([0-9]+\.[0-9]+)\;/){ $8=substr($0, RSTART, RLENGTH-1) }1' | awk '{ sub("QD=","",$8); print }' | awk '$8>7'
# Or
bcftools query -f '[%ALT\t%QD\n]' files/sample.vcf.gz | awk '$2>7'

#10. How many contigs are referred to in the file. Check the header section
zcat files/sample.vcf.gz | grep -c '^##contig'

#11. Comment on the eighth and ninth columns of the file
zcat files/sample.vcf.gz |  sed -n '/^[#?][^#]/, $p' | cut -f 8,9

#12. Extract data on the read depth of called variants for sample SRR13107018
bcftools query -f '%ALT\t%DP\n' -s SRR13107018 files/sample.vcf.gz

#13. Extract data on the allele frequency of alternate alleles. Combine this data with the chromosome and position of the alternate allele
bcftools query -f '%CHROM\t%POS\t%ALT\t%AF\n' files/sample.vcf.gz

## SAM file manipulation

#1. Describe the format of the file and the data stored
## Creating required bam files
samtools view -bS files/sample.sam > files/sample.bam
samtools sort files/sample.bam > files/sample_sorted.bam
samtools index files/sample_sorted.bam

#2. What does the header section of the file contain
samtools view -H files/sample.sam

#3. How many samples are in the file
grep -c SM files/sample.sam

#4. How many alignments are in the file
awk ' $1 !~ /@/ {print $1}' files/sample.sam | wc -l

#5. Get summary statistics for the alignments in the file
samtools view -bS files/sample.sam | samtools stats > out_dir/samstats.txt

#6. Count the number of fields in the file
grep -v "^@" files/sample.sam | awk '{print NF}'| sort -nu

#7. Print all lines in the file that have @SQ and sequence name tag beginning with NT_
grep "@SQ.*NT_" files/sample.sam

#8. Print all lines in the file that have @RG and LB tag beginning with Solexa
grep "@RG.*LB:Solexa" files/sample.sam

#9. Extract primarily aligned sequences and save them in another file
samtools view -F 4 files/sample.sam  > aligned_reads

#10. Extract alignments that map to chromosomes 1 and 3. Save the output in BAM format
samtools view -b files/sample_sorted.bam "1:3" > out_dir/align_1_3.bam

#11. How would you obtain unmapped reads from the file
samtools view -f 4 files/sample.bam

#12. How many reads are aligned to chromosome 4
grep -c "^4\t" sample.sam

#13. Comment of the second and sixth column of the file
## Comments in README.md file.

#14. Extract all optional fields of the file and save them in “optional_fields.txt
awk '{for(i=11;i<=NF;i++) print $i}' files/sample.sam > out_dir/optional_fields.txt
