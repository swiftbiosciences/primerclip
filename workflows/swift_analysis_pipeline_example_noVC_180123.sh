#!/bin/bash

####################### EXAMPLE ANALYSIS WORKFLOW ############################
# NOTES: Basic steps/workflow to analyze fastq files:
#   1. adapter trimming
#   2. alignment
#   3. primerclip
#   4. sam sort and index

# * PICARD used as default. SAMTOOLS example also available. Comment out PICARD
# and use SAMTOOLS.

# check that a masterfile as been given as command line argument
# if not, print the help message to stdout
missing=$'Swift Biosciences Example Analysis Script
        Missing: MASTERFILE
                 PLEASE INCLUDE MASTERFILE AS ARGUMENT
        Useage: ./swift_analysis_pipeline_example_withVC_180118.sh MASTERFILE

        Description: Start with fastq files and get coverage plus
                     on-target metrics in addition to vcf files
                     from two variant callers.'

if [ $# == 0 ]; then
    echo "$missing"
    exit 1
fi

echo "Swift Biosciences Analysis Workflow without Metrics + Variant Calling"
# reference genome - human
ref=/seq/refgenomes/Homo_sapiens_assembly19broad.fasta
picard='java -jar /tools/picard/dist/picard.jar'
# panel masterfile
master="$1"
echo "accepted masterfile"

# ANALYSIS WORKFLOW
# Repeat steps for each sample
for f in *_R1_001.fastq.gz
do

# assign fastq R1 and R2 files and prefix of sample name to variables
fq1=$f
fq2=${fq1%%_R1_*}_R2_001.fastq.gz
prefix=${fq1%%_R1_001.fastq.gz}

# ADAPTER TRIMMING
# trimmomatic36
trimdir=/tools/trimmomatic36

echo "trimming Illumina adapters"
java -Xmx24g -Xms16g -jar ${trimdir}/trimmomatic-0.36.jar PE \
    -threads 12 -trimlog ${prefix}_trimmatic_trimlog.log \
    $fq1 $fq2 ${prefix}_R1_adapter_trimmed.fq.gz ${prefix}_unpaired_R1.fq.gz \
    ${prefix}_R2_adapter_trimmed.fq.gz ${prefix}_unpaired_R2.fq.gz \
    ILLUMINACLIP:${trimdir}/adapters/TruSeq3-PE.fa:2:30:10 \
    MINLEN:30 2> ${prefix}_adapter_trimming.log

fqt1=${prefix}_R1_adapter_trimmed.fq.gz
fqt2=${prefix}_R2_adapter_trimmed.fq.gz

# READ ALIGMENTS
# bwa
echo "aligning with bwa b37"
bwa mem $ref $fqt1 $fqt2 -U 17 -M -t 32 > ${prefix}_non_primer_trimmed.sam \
    2> ${prefix}_bwa_alignment.log

# PRIMER TRIMMING
# primerclip
echo "trimming primers"
primerclip-lowmem-v171205 ${prefix}_non_primer_trimmed.sam "${master}" \
    ${prefix}_primer_trimmed.sam

# SAM TO BAM FILE OPTIONS
# TWO OPTIONS to sort and index sam to bam files: 
#   1. with picard or 2. without picard
#   **Both options shown for primer and non primer trimmed sam files**

# 1. WITH PICARD: sorted and indexed sam to bam file
# non primer trimmed sam file to bam file
echo "sorting and indexing non primer trimmed sam file"
$picard SortSam I=${prefix}_non_primer_trimmed.sam O=${prefix}_non_primer_trimmed.bam \
    CREATE_INDEX=true SORT_ORDER=coordinate \
    VALIDATION_STRINGENCY=LENIENT 2> ${prefix}_make_non_primer_trimmed_bam.log

samtools index ${prefix}_non_primer_trimmed.bam 2> ${prefix}_index_nonptrim.log

# primer trimmed sam file to bam file
echo "sorting and indexing primer trimmed sam file"
$picard SortSam I=${prefix}_primer_trimmed.sam O=${prefix}_primer_trimmed.bam \
    CREATE_INDEX=true SORT_ORDER=coordinate \
    VALIDATION_STRINGENCY=LENIENT 2> ${prefix}_make_primer_trimmed_bam.log

samtools index ${prefix}_primer_trimmed.bam 2> ${prefix}_index_ptrim.log

# 2. WITHOUT PICARD: unsorted sam to bam file
# non primer trimmed sam file to bam file
#echo "sorting and indexing non primer trimmed sam file to bam file"
#samtools view -S -b ${prefix}_non_primer_trimmed.sam \
#    > ${prefix}_non_primer_trimmed.bam

#samtools sort ${prefix}_non_primer_trimmed.bam -o ${prefix}_non_primer_trimmed_sorted.bam
#samtools index ${prefix}_non_primer_trimmed_sorted.bam
# primer trimmed sam file to bam file
#echo "sorting and indexing primer trimmed sam file to bam file"
#samtools view -b -S ${prefix}_primer_trimmed.sam \
#    > ${prefix}_primer_trimmed.bam

#samtools sort ${prefix}_primer_trimmed.bam -o ${prefix}_primer_trimmed_sorted.bam
#samtools index ${prefix}_primer_trimmed_sorted.bam

echo "completed"
# check markdown for code to convert from bam to sam
rm *.sam
done
