# vcf_sam_manipulation_assn

## Bash: Manipulation of the VCF and SAM files

### By : KAKANDE Paul

# VCF files.

## 1-2. Structure of VCF file and contents of the header.

> The vcf file is made of _two_ main parts:
> - Header
> - Variant Call Records
>
> ### Header
>
> The header lines start with **##**.
> Its first line indicates the version of the vcf format here given as VCFv4.2.
> Next are the filter lines indicating the filters applied to the data.
>> For this sample file, the filter used is `LowQual` 
> Next are the format lines `##FORMAT`
> These lines are formated with descriptors **ID** for the short forms used in variant call records, **Type** indicating the type of data that is whether Integer, String among others, and **Description** which provides more information about the annotation in question.
>> These lines provide explanations for the different annotations used in the Variant Call Records section for example:
>> - AD : Allelic depth fo the ref and alt alleles in the order listed
>> - DP : Appoximate read depth (reads with MQ=255 or with bad mates are filtered)"
>> - among others.
> The next lines contain the GATKCommnadLine `##GATKCommandLIne`
> These lines contain the parameters that were used by the tool that generated the file.
>> In this sample file, the tool is `GenomicsDBImport`, `GenotypeGVCFs`, and `HaplotypeCaller`
> The final lines are the Contig lines and Reference.
> These contain the contig names, lengths and which reference assembly used.

> ### Variant Call Records
>
> This section is tab delimited.
> The first line is descriptive of the columns with the first **8** columns representing the properties observed at the level of the variant (or invariant) site. The first **7** are required by VCF format but may be empty filled with a **.** as a placeholder. These include 
> - **CHROM and POS** : These annotate for the contig and genomic coordinates on which the variant occurs.
> - **ID** : Optional identifier for the variant.
> - **REF and ALT** : Reference and alternative allele(s) observed.
> - **QUAL** : The Phred-scaled probability that a REF/ALT polymorphism exists at this site given sequencing data. Because the Phred scale is -10 * log(1-p), a value of 10 indicates a 1 in 10 chance of error, while a 100 indicates a 1 in 10^10 chance .
> - **FILTER** : This contains the name(s) of the filters that the variant fails to pass or the value PASS if it passed all filters. If the FILTER contains `.`, then no filtering has been applied to the records.
>
> The eighth column, the **INFO** column is not a requirement for VCF.
> The information contained in here is represented as tag-value pairs with the tag and value separated by an equal sign.
> Their description is provided for in the header section `FORMAT` part.

> The next columns are sample-level annotations.
> `FORMAT`  and then followed by the different samples.
 
> With the chromosomes, positions, ID, Reference and Alternate alleles, qual, filter, information on allele counts, frequency, depth, their format, and 6 columns of 6 samples. In total, the vcf file has 15 columns. 

Commands used:
`zcat files/sample.vcf.gz` to view the whole file
`zcat files/sample.vcf.gz |  grep '^##'` to view the header in isolation.

Or 
`bcftools view files/sample.vcf.gz` to view the whole file
`bcftools view -h files/sample.vcf.gz` to view the header in isolation.

## 3. Number of Samples in file.
Command used:
`zcat files/sample.vcf.gz | grep -m1 "^#CHROM" | cut -f 10- | tr "\t" "\n"  | wc -l`

Or 
`bcftools query -l files/sample.vcf.gz | wc -l`

The output is **6** indicating six samples are in the file

## 4. Number of Variants in the file.
Command used:
`zcat files/sample.vcf.gz |  grep -c '^[^#]'`

Or 
`bcftools query -f '%ALT\n' files/sample.vcf.gz | wc -l`

The output is **398246**.

## 5. Extraction of  the chromosome, position, QualByDepth and RMSMappingQuality fields and save the output to a tab-delimited file.
Commmand used:
`zcat files/sample.vcf.gz |  grep '^[^#]' | cut -f 1,2,6,8 | awk 'match($0, /MQ=([0-9]+\.[0-9]+)\;/){ $4=substr($0, RSTART, RLENGTH-1) }1' | awk '{ sub("MQ=","",$4); print }' OFS='\t' > out_dir/answer5.tsv`

Or 
`bcftools query -f '%CHROM\t%POS\t%QUAL\t%MQ\n' files/sample.vcf.gz > out_dir/answers5.tsv`

## 6. Data belonging to chromosomes 2,4 and MT
Command used:
`zcat files/sample.vcf.gz | awk -F "\t" '$1~/^(2|4|MT)$/{print $0}'`

## 7. Variants that do not belong to chr20:1-30000000
Command used:
`zcat files/sample.vcf.gz |  grep '^[^#]' | awk '!($1 == 20 && $2 >= 1 && $2 <= 30000000)'`

## 8. Variants that belong to SRR13107019
Command used:
`bcftools query -f '%ALT\n' -s SRR13107019 files/sample.vcf.gz`

## 9. Variants with a QuadByDepth above 7
Command used:
`zcat files/sample.vcf.gz |  grep '^[^#]' | awk 'match($0, /QD=([0-9]+\.[0-9]+)\;/){ $8=substr($0, RSTART, RLENGTH-1) }1' | awk '{ sub("QD=","",$8); print }' | awk '$8>7' | awk 'NF > 0'|wc -l`

Or 
`bcftools query -f '%ALT\t%QD\n' files/sample.vcf.gz | awk '$2>7'`

## 10. Number of contigs
Command used:
`zcat files/sample.vcf.gz | grep -c '^##contig'`

The number of contigs obtained was **2211**

## 11. Comment on eight and ninth columns
These were viewed using `zcat files/sample.vcf.gz |  sed -n '/^[#?][^#]/, $p' | cut -f 8,9`

>The eighth column; the INFO column stores information about the variant itself. It contains information such as the variant's quality (QUAL), read depth (DP), and allele frequency (AF). The format of this column is not fixed, that is different tools may add more fields such as below:
> - AC=2;AF=1.00;AN=2;DP=2;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=21.13;SOR=0.693
> - AC=8;AF=0.667;AN=12;BaseQRankSum=1.25;DP=145;ExcessHet=1.3830;FS=0.000;MLEAC=8;MLEAF=0.667;MQ=60.00;MQRankSum=0.00;QD=28.73;ReadPosRankSum=1.34;SOR=0.928 

>The ninth column; FORMAT column contains information about the genotype of the sample(s) at a specific variant location. It contains information such as the genotype (GT), read depth (DP), and allelic depth (AD) for each sample. The format of this column is fixed and defined in the VCF header such as:
> - GT:AD:DP:GQ:PL 

## 12. Read depth variants for sample SRR13107018
Command used:
`bcftools query -f '%ALT\t%DP\n' -s SRR13107018 files/sample.vcf.gz`

## 13. Allele frequency of alternate alleles extract combined with chromosome postion and alternate allele
Command used:
`bcftools query -f '%CHROM\t%POS\t%ALT\t%AF\n' files/sample.vcf.gz`


# Manipulating SAM files
## 1. Describe the format of the file and the data stored
The Sequence Alignment (SAM) file is a tab delimited text file that stores information about the alignment of reads in a FASTQ file to a reference genome or transcriptome.
A SAM file format is consists of:
- Header that contains sequence dictionary, read group definitions among others.
- Records containing the structured read information one line per record


## 2. What does the header section of the file contain
Command used:
`samtools view -H files/sample.sam`.

>The header section contains; 
>- HD- general information about the alignment file with fields: 
>-- GO- group ordering 
>-- SO- sort ordering of the alignments in the file
>- VN- the version number of the header 
>- SQ- information about  the reference sequence used for the alignment 
>- SN- name of the reference sequence
>- LN- length of the reference sequence.

>The next part of the header contains the :
>- RG- read groups from the same sample or sequencing run
>- ID- identifies the read group
>-  PL-sequencing platform used to generate read groups
>- PU-specific sequencing platform unit from which the reads were generated
>- LB- stores the name of the library from which the reads were generated 
>- SM- name of the sample the reads come from
>- CN- sequencing center where the reads were generated. 

>Other information contained in the header is BI which is barcode indexing during multiplexing and PE for paired ended reads.

>The bottom of the header contains
>- PG- the program that generated the alignment
>- PN-program name and VN- the program version
>- PP- the previous program used before the alignment
>- ID-the command line argument used to run the program that generated the alignment. 

## 3. How many samples are in the file
Command used:
`grep -c SM files/sample.sam`

The output is **249** refering to samples.

## 4. How many alignments are in the file
Command used:
`awk ' $1 !~ /@/ {print $1}' files/sample.sam | wc -l`

The output is **361428* refering to alignments.

## 5. Get summary statistics for the alignments in the file
Command used:
`samtools view -bS files/sample.sam | samtools stats > out_dir/samstats.txt`

The sam file is converted in a Binary Alignment Map (BAM) using the `samtools view -bS files/sample.sam` then the output is piped to `samtools stats` whose output is stored in file **samstats.txt** in the out_dir directory.

## 6. Count the number of fields in the file
Command used:
`grep -v "^@" files/sample.sam | awk '{print NF}'| sort -nu`

The first part removes the header and awk is used to count the number of fields.
Three values are returned that is **13,17 and 18** as some lines have more or less fields than others.

## 7. Print all lines in the file that have @SQ and sequence name tag beginning with NT_
Command used:
`grep "@SQ.*NT_" files/sample.sam`

## 8. Print all lines in the file that have @RG and LB tag beginning with Solexa
Command used:
`grep "@RG.*LB:Solexa" files/sample.sam`

## 9. Extract primarily aligned sequences and save them in another file

# 10. Extract alignments that map to chromosomes 1 and 3. Save the output in BAM format

## 11. How would you obtain unmapped reads from the file
## 12. How many reads are aligned to chromosome 4
## 13. Comment of the second and sixth column of the file
## 14. Extract all optional fields of the file and save them in “optional_fields.txt
