#!/bin/bash

#Author: Osman Merdan
#Date Last Modified: 04/12/2022

# Stop on any error 
# -e: Exit immediately if a command exits with a non-zero status.
# -u  Treat unset variables as an error when substituting.
# -x  Print commands and their arguments as they are executed.
set -ue

#####################################################################################################################
#                                             1-VARIANT CALLING                                                     #
#####################################################################################################################

# List of SRA accession numbers are obtained from:
# https://www.ncbi.nlm.nih.gov/Traces/study/?acc=ERP132402&o=acc_s%3Aa

# ...........................#
# 1-Downloading SRA sequences
# ...........................# 
mkdir -p "raw-data"
parallel --bar --verbose "fastq-dump --split-files --outdir ./raw-data {}" < ./SRR_Acc_List.txt

# ..........................................#
# 2-Trimming and removing low quality reads 
# ..........................................#
# From the Trimmomatic Documentation ;
# 1. ILLUMINACLIP
#  ./TruSeq3-PE.fa: Adapter sequences.
#  Required adapter sequences can be found at:
#  https://github.com/usadellab/Trimmomatic/tree/main/adapters
#  :2 = seedMismatches: 
#  specifies the maximum mismatch count which will still allow a full match to be performed
#  Initially Trimmomatic will look for seed matches (16 bases) allowing maximally 2 mismatches.
#  :30 = palindromeClipThreshold:
#  how accurate the match between the two 'adapter ligated' reads must be for PE palindrome read alignment. 
#  :10 = simpleClipThreshold:
#  specifies how accurate the match between any adapter etc. sequence must be against a read.
#  :2 = minAdapterLength: 
#  In addition to the alignment score, palindrome mode can verify that a minimum length of adapter has been detected. 
#  If unspecified, this defaults to 8 bases, for historical reasons.
#  However, since palindrome mode has a very low false positive rate -->
#  this can be safely reduced, even down to 1, to allow shorter adapter fragments to be removed.
#  True: Always keep both reads. : The reverse read will also be retained in palindrome mode.
#  Dropping reverse reads can cause problems in downstream analysis. 
# 2. SLIDINGWINDOW
#  :7 = window size
#  :25 = required quality
mkdir -p "paired"
mkdir -p "unpaired"
while read -r i;
do
trimmomatic PE\
 ./raw-data/"${i}"_1.fastq\
 ./raw-data/"${i}"_2.fastq\
 ./paired/"${i}"_paired_1.fastq\
 ./unpaired/"${i}"_unpaired_1.fastq\
 ./paired/"${i}"_paired_2.fastq\
 ./unpaired/"${i}"_unpaired_2.fastq\
 ILLUMINACLIP:./TruSeq3-PE.fa:2:30:10:2:True SLIDINGWINDOW:7:25 LEADING:20 MINLEN:200
done < ./SRR_Acc_List.txt

# Zip raw data folder and its contents. 
# From now on every analysis will be performed on trimmed paired end dataset
gzip ./raw-data/*.fastq


# ............................................#
# 3-Quality Assesment of Trimmed Raw Sequences 
# ............................................#
mkdir -p ./qc
parallel --bar --verbose "fastqc {} -o ./qc" ::: ./paired/*.fastq
multiqc ./qc -o ./qc

# ............................................#
# 4-Preparing Files For Alignment
# ............................................#
mkdir -p ref 
mkdir -p alingments
bio fetch AY446894.2 > ./ref/cmv_merlin.gb # Genbank record 
bio fasta ./ref/cmv_merlin.gb > ./ref/cmv_merlin.fa # Converting fasta
samtools faidx ./ref/cmv_merlin.fa # Indexing genome for alignment
bwa index ./ref/cmv_merlin.fa # Indexing for IGV browser visualization

# ............................................#
# 5-Alignment (BWA-MEM)
# ............................................#
while read -r i; do
bwa mem -t 6\
 ./ref/cmv_merlin.fa\
 ./paired/"${i}"_paired_1.fastq\
 ./paired/"${i}"_paired_2.fastq > ./alignments/"${i}".sam
done < ./SRR_Acc_List.txt

# ............................................#
# 6-Post Proccessing Alignment Files 
# ............................................#
# 1. fixmate -r
# Removes unmapped reads and secondary alignments
# 2. sort
# Sorts alignments according to chr coordinates
# 3. wiev for filtering
#  -q 30: min-map quality 30
#  -F 2048: exclude 2048 flag (supplementary alignment)
#  -f 2: select only proper pairs (flag 2) 
#  -h: with header
mkdir -p bam
parallel --verbose --bar\
 "samtools fixmate -O bam -r {} -\
 | samtools sort -\
 | samtools view -q 30 -F 2048 -f 2 -b -h -O bam -o ./bam/{/.}.bam"\
 ::: ./alignments/*
 
# ............................................#
# 7-Duplicate removal
# ............................................#
mkdir -p deduplicated-bam
parallel --verbose --bar\
 "picard MarkDuplicates --INPUT {}\
 --OUTPUT ./deduplicated-bam/{/.}.bam\
 --REMOVE_DUPLICATES true\
 --ADD_PG_TAG_TO_READS false\
 --METRICS_FILE ./deduplicated-bam/{/.}-metrics.txt\
 --REMOVE_SEQUENCING_DUPLICATES true"\
 ::: ./bam/*.bam

# Indexing bam files for variant calling
parallel --bar --verbose\
 "samtools index {}" ::: ./deduplicated-bam/*.bam

# ............................................#
# 8-Alignment Stats 
# ............................................#
mkdir ./alignment-stats
parallel --verbose --bar\
 "samtools flagstat {} > ./alignment-stats/{/.}-flagstat.txt"\
 ::: ./deduplicated-bam/*.bam

# Summarizing alignment stats with MultiQC
multiqc ./alignment-stats -o ./alignment-stats

# Coverage stats
# Pipe friendly tabular format general coverage stats per sample
mkdir -p ./coverage
printf\
 "name\tstart\tend\tnumreads\tcovbasses\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq\tsampleID\n"\
 > ./coverage/coverage-stats.txt

while read -r name;
do
samtools coverage --no-header ./deduplicated-bam/"${name}".bam |\
 awk -F '\t'\
 -v name="$name"\
 '{print $0 "\t" name }'\
 >> ./coverage/coverage-stats.txt
done < ./SRR_Acc_List.txt

# Depth stats
parallel --verbose --bar\
 "samtools depth -a -o ./coverage/{/.}.txt {}"\
 ::: ./deduplicated-bam/*.bam

# Concatenating depth of coverage-position files
#  -F '\t': input file is tab delim 
#  -v i="$i": copies shell variable in awk program
#  $i = <<SRR Accession Number>> = i variable in awk program
#  print $2 "\t" $3 "\t"i --> prints tab seperated column 2, 3 and SRR Accession 
#  position + depth + <<SRR Accession Number>>  
while read -r i; 
do
      awk -F '\t'\
      -v i="$i"\
      '{print $2 "\t" $3 "\t"i}'\
      ./coverage/"$i".txt\
      >> ./coverage/piled-depth.txt 
done < ./SRR_Acc_List.txt

# ............................................#
# 9- Getting Rid of Low-Cov/Depth Alignments 
# ............................................#
# In order to achieve relaible variant calling results: 
# Mean depth requirement was determined as and 
# Overall genome coverage was determined as 90%
# Those variables were stored in coverage-stats.txt
# Using awk command to find sampleIDs which have :
# coverage>90 or mean-depth<20
# Printing sampleIDs {column 10} to a new txt file
mkdir -p ./deduplicated-bam/low-meandepth-samples
awk -F "\t"\
 '($6<94 || $7<20) && NR>1 {print $10}'\
 ./coverage/coverage-stats.txt\
 > ./deduplicated-bam/low-meandepth-samples/low-meandepth-samples.txt

# Using the low-meandepth sampleIDs,
# move low-meandepth bam files to a new directory. 
while read -r i;
do
mv ./deduplicated-bam/"${i}".bam\
 ./deduplicated-bam/"${i}".bam.bai\
 ./deduplicated-bam/"${i}"-metrics.txt\
 ./deduplicated-bam/low-meandepth-samples
done < ./deduplicated-bam/low-meandepth-samples/low-meandepth-samples.txt

# Storing high quality sample names to new txt file for further processing.
# Sorted files is necessary to perform this opeation.
# comm function -->
#  -1: Suppress printing of column 1, lines only in file1.
#  -3: Suppress printing of column 3, lines common to both.
#  -i: Case insensitive comparison of lines.
comm -13i\
 ./deduplicated-bam/low-meandepth-samples/low-meandepth-samples.txt\
 ./SRR_Acc_List.txt\
 >./SRR_further_processed.txt

# ............................................#
# 10-Variant calling and filtering 
# ............................................#
# 1. bcftools mpileup
#  -O u: output type uncompressed bcf for memory efficiency
#  -d 2000: At a position, read maximally 2000 reads per input file.
#  -a: Comma-separated list of FORMAT and INFO tags to output.
#  -f: The faidx-indexed reference file in the FASTA format.
#    'INFO/AD': Total allelic depth  
#    'INFO/ADF': Total allelic depths on the forward strand
#    'INFO/ADR': Total allelic depths on the reverse strand
#    'FORMAT/DP':Number of high-quality bases
# 2. bcftools call
#  --ploidy 1: predefined ploidy
#  --variants-only: output variant sites only 
#  --keep-alts: output all alternate alleles present in the alignments
#  -m: Alternative model for multiallelic and rare-variant calling
#  --output-type v: uncompressed vcf
# 3. bcftools norm: 
# Left-align and normalize indels; check if REF alleles match the reference;
# split multiallelic sites into multiple rows; recover multiallelics from multiple rows.
#  -m: Split multiallelics (-) or join biallelics (+) -|+[snps|indels|both|any]
#  -f: Refrence file 
# 4. bcftools filter
#  -i: include
#     (AD[1]/(FMT/DP)>0.05): Alternate alele/number of high quality bases ratio
#     QUAL>30: Variant quality>30
#     FS<10 : Fisher Strand (phred scales p value for strans bias)
#             phread score 10 --> p value = 0.1 
#             p-value of 1 (FS=0) meaning there is a 100% chance of there being no bias.
#             More detail here : 
#             https://gatk.broadinstitute.org/hc/en-us/articles/360035532152-Fisher-s-Exact-Test
#     AD[1]>10: At least 10 observations for alternete allele
#     FORMAT/DP>20: Number of high-quality bases at least 20 at a given position
mkdir -p vcf
parallel --bar --verbose\
 "bcftools mpileup\
 -O u\
 -d 2000\
 -f ./ref/cmv_merlin.fa\
 -a 'INFO/AD,INFO/ADF,INFO/ADR,FORMAT/DP' {}\
 | bcftools call --ploidy 1 --variants-only --keep-alts -m --output-type v\
 | bcftools norm -f ./ref/cmv_merlin.fa -m -any\
 | bcftools filter -i '(AD[1]/(FORMAT/DP)>0.05) && QUAL>30 && FS<10 && AD[1]>10 && FORMAT/DP>20'\
 > ./vcf/{/.}.vcf"\
 ::: ./deduplicated-bam/*.bam

# ............................................#
# 11-Variant annotation 
# ............................................#
# Building SnpEff Database 
# Note: snpEff package has a built in script to download genomes 
# In conda enviroment the script is stored in below path
# ~/miniconda3/envs/<envname>/share/snpeff-<version>/scripts/buildDbNcbi.sh
# The script need only one argument : GenBank accession number
# Don't forget to go in snpeff folder before building database!!
# cd  ~/miniconda3/envs/<envname>/share/snpeff-<version>/ then run the script.
# Note: SnpEff build returns error: 
# Transcript 'HHV5wtgr002' is already in Gene 'HHV5wtgr002'
# Workaround for this issue is remove duplicate mrna entries in genbank file.
# And running the build function manually.

# Annotating variants and saving summary htmls 
mkdir -p ./annotated-vcf
mkdir -p ./html-vcf-stats
parallel --bar --verbose\
 "snpEff eff\
 -no-downstream\
 -no-intergenic\
 -no-intron\
 -no-upstream\
 -no-utr\
 -htmlStats ./html-vcf-stats/{/.}.html\
 AY446894.2 {} > ./annotated-vcf/{/.}.vcf"\
 ::: ./vcf/*.vcf

# ............................................#
# 12- Annotated variant extraction
# ............................................#
# Extracting fields using SnpSift into tabular format for downstream analysis
# To extract fields one per line useing the script from:
# https://github.com/pcingola/SnpEff/blob/master/scripts/vcfEffOnePerLine.pl 
# Further documentation here https://pcingola.github.io/SnpEff/ss_extractfields/
wget\
 https://raw.githubusercontent.com/pcingola/SnpEff/master/scripts/vcfEffOnePerLine.pl \
 -P ./
mkdir -p extracted-vcf
parallel --bar "cat {}\
 |perl  ./vcfEffOnePerLine.pl\
 |SnpSift extractFields - -s ',' -e '.' 'CHROM'\
 'POS' 'ID' 'REF' 'ALT' 'FILTER' 'ANN[*].EFFECT'\
 'ANN[*].IMPACT' 'ANN[*].GENE' 'ANN[*].GENEID'\
 'ANN[*].FEATURE' 'ANN[*].FEATUREID'\
 'ANN[*].HGVS_C' 'ANN[*].HGVS_P' > ./extracted-vcf/{/.}" ::: ./annotated-vcf/*.vcf

# ............................................#
# 13- Concatenate extracted variants
# ............................................#
# Concatenate all extracted vcf files to one to work on later
# Note that non coding variants were not annotated. (EFF=.)
# $7!="." --> removes records without annotation
# NR!=1 --> removes first line of the file (column names)
# sub(/.*\//, "", FILENAME) --> removes file path from FILENAME
# print $0 "\t"FILENAME --> puts tab+FILENAME at the end of each row.
for file in ./extracted-vcf/*
do
     awk -F '\t'\
     '($7!=".") && (NR!=1)\
     {sub(/.*\//, "", FILENAME); print $0 "\t"FILENAME}'\
     "$file"\
      >> ./extracted-vcf/merged-variants.txt
done




# -------------------------------------------------------------------------------------------------------------------




#####################################################################################################################
#                                                  2-GENOTYPING                                                     #
#####################################################################################################################

# ............................................#
# 1- Motif - Read Matching
# ............................................#
# Following similar steps described in below papers 
# https://doi.org/10.3390/v14050855
# https://doi.org/10.1093/infdis/jiz208
mkdir -p ./genotyping
# Downloading miRNA_Search.c from Centre for Virus Research Github page
# MIRNA_SEARCH finds sequence patterns in fasta files 
wget\
 https://github.com/centre-for-virus-research/VATK/blob/master/GenotypingTools/miRNA_Search.c\
 -P ./genotyping
wget\
 https://github.com/centre-for-virus-research/VATK/blob/master/GenotypingTools/miRNA_Search.h\
 -P ./genotyping
# Install MIRNA-Search 
gcc ./genotyping/miRNA_Search.c -o ./genotyping/MIRNA_SEARCH
# Downloading motifs from Centre for Virus Research Github page
wget\
 https://github.com/centre-for-virus-research/VATK/raw/master/HCMV_pipeline/sig-12Genes_13Sep17.fa\
 -P ./genotyping



















#################################
#  VCF file to Consensus Fasta  #
#################################

# bcftools consensus creates consensus fasta files for samples
# seqkit replace replaces fasta headers with corresponding sample name
# Note : In this approach indels were taken into account
# If the main interest is SNP's only -->
# Follow theese steps :
# 1-Filter out indels using vcftools
# 2-Create consensus and cocatenate consensus files 
# 3-Concatenated consensus file is equivalent to MSA
# 4-Find SNP sites which is equivalent to MSA
# 5-Create RAxML tree
mkdir -p consensus
bgzip ./vcf/*.vcf
bcftools index ./vcf/*.vcf.gz
while read -r name;
do
bcftools consensus ./vcf/"${name}".vcf.gz\
 --fasta-ref /Users/osmanmerdan/cmv-project/ref/cmv_merlin.fa\
 seqkit replace --pattern '.+'\
 --replacement "${name}"\
 > ./consensus/"${name}".fa
done < ./SRR_further_processed.txt
# Concatenate consensus files and refrence file into one big fasta file 
cat ./consensus/*.fa ./ref/cmv_merlin.fa > one-big.fa


####################
#  MSA with MAFFT  #
####################

# MAFFT uses lots of computing power in accuracy modes. 
# Thats why method was set to auto
# --reorder : sort fasta outout files
# --treeout : creates phylogenetic tree 
 mafft --auto --reorder --treeout one-big.fa > msa.fa
# Note : This program will spit out one-big.fa.tree 
# That tree file can be visualized using R ape package.

#############################################
# Phylogenetic Tree and SNP Distance Matrix #
#############################################

# snp-sites finds snp sites from a multi FASTA alignment file.
# -p --> output a phylip file (for RAxML and other tree building applications.)
# -o STR specify an output filename
snp-sites -p -o "snp-sites.phy" msa.fa

# RAxML (Randomized Axelerated Maximum Likelihood) -->
# https://cme.h-its.org/exelixis/web/software/raxml/hands_on.html
# PThreads --> Parallel, if you use the PThreads version they will run faster!
# -f??a --> rapid??Bootstrap??analysis??and??search??for??bestscoring??ML??tree??in??one??program??run
# -m GTRGAMMA --> ??m??substitutionModel,  ML search on DNA data
# As starting trees RAxML uses randomized stepwise addition parsimony trees -->
## it will not generate the same starting tree every time.
## You can either force it by providing a fixed random number seed --> -p 12345
# -x 12345 rapid bootstrap random number seed
# -T -->??numberOfThreads
# -#100 --> 100 rapid Bootstrap searches
# This RAxML will do 100 rapid Bootstrap searches,
# 20 ML searches and return -->
# the best ML tree with support values to you via one single program call.
# -s sequenceFileName
# -n outputFileName
raxmlHPC-PTHREADS -T 6 -f a -m GTRGAMMA -p 12345 -x 12345\
 -#100\
 -s snp-sites.phy\
 -n ml-tree

# snp-dists --> Pairwise SNP distance matrix from a FASTA alignment
# -j --> number of threads 
# -a --> count all differences 
# -b --> blank top left corner cell
snp-dists -j 6 -a -b msa.fa > ./matrix.tsv 



######
#Compress Remainin Files 
######
#tar -czvf name-of-archive.tar.gz /path/to/directory-or-file
#c: Create an archive.
#-z: Compress the archive with gzip.
#-v: Display progress in the terminal while creating the archive, also known as ???verbose??? mode. The v is always optional in these commands, but it???s helpful.
#-f: Allows you to specify the filename of the archive.