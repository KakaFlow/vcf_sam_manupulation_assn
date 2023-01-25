# vcf_sam_manupulation_assn

## Bash: Manupulation of the VCF and SAM files

## By : KAKANDE Paul

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
`bcftools view files/sample.vcf.gz` to view the whole file
`bcftools view -h files/sample.vcf.gz` to view the header in isolation.

## 3. Number of Samples in file.
Command used:
`bcftools query -l files/sample.vcf.gz | wc -l`
The output is **6** indicating six samples are in the file

## 4. Number of Variants in the file.
Command used:
`bcftools query -f '%ALT\n' files/sample.vcf.gz | wc -l`
The output is **398246**.

## 5. Extraction of  the chromosome, position, QualByDepth and RMSMappingQuality fields and save the output to a tab-delimited file.
Commmand used:
`bcftools query -f '%CHROM\t%POS\t%QUAL\t%MQ\n' files/sample.vcf.gz > out_dir/answers5.tsv`

## 6. Data belonging to chromosomes 2,4 and MT
Command used:


