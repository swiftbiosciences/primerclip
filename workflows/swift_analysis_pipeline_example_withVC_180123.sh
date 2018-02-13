#!/bin/bash
####################### EXAMPLE ANALYSIS WORKFLOW ############################
# NOTES: Basic steps/workflow to analyze fastq files:
#   1. adapter trimming
#   2. alignment
#   3. primerclip
#   4. sam sort and index
#   5. calculate metrics: coverage + on-target
#   7. variant calling with GATK + LoFreq

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

echo "Swift Biosciences Analysis Workflow with Metrics + Variant Calling"
# reference genome - human
ref=/seq/refgenomes/Homo_sapiens_assembly19broad.fasta
# tools' variables
picard='java -jar /tools/picard/dist/picard.jar'
gatk='java -jar /tools/GenomeAnalysisTK.jar'

# panel masterfile
echo "accepted masterfile"
master="$1"

# panel bedfiles from masterfile
# used to calculate coverage and on target metrics
echo "making bedfiles from masterfile"
# non merged target bedfile
awk '{print $1,$2,$3,$4}' OFS="\t" "${master}" |
    sort -k1,1n -k2,2n > nonmerged_targets.bed
# merged targets bedfile
bedtools merge -nms -i nonmerged_targets.bed | sed 's/;.*//' | 
    sort -k1,1n -k2,2n > merged_targets.bed
# merged targets bed in 5 column format for bedtools compatibility
awk '{print $1,$2,$3,"+",$4}' OFS="\t" merged_targets.bed \
    > merged_targets_5col.bed

# ANALYSIS WORKFLOW
# Repeat analysis steps for each sample
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

# SAM TO BAM FILE
# picard tools
# primer trimmed sam file to primer trimmed bam file
echo "sorting and adding read groups"
# AddOrReplaceReadGroups needed to assign annotation to SAM file for
# GATK to function correctly. Can change the tags as needed by user.
# ie. hiseq instead of miseq etc
$picard AddOrReplaceReadGroups \
    I=${prefix}_primer_trimmed.sam O=${prefix}_primer_trimmed.bam \
    SO=coordinate RGID=snpID LB=swift SM=${prefix} PL=illumina PU=miseq \
    2> ${prefix}_add_readgroups.log
# indexing primer trimmed bam file
samtools index ${prefix}_primer_trimmed.bam 2> ${prefix}_idx_ptrimed_bam.log

# non primer trimmed sam file to non primer trimmed bam file for debugging
$picard SortSam I=${prefix}_non_primer_trimmed.sam O=${prefix}_non_primer_trimmed.bam \
    CREATE_INDEX=true SORT_ORDER=coordinate \
    VALIDATION_STRINGENCY=LENIENT 2> ${prefix}_make_non_primer_trimmed_bam.log
# indexing non primer trimmed bam file
samtools index ${prefix}_non_primer_trimmed.bam 2> ${prefix}_idx_nonptrimmed_bam.log

# CALCULATING COVERAGE METRICS AND ON TARGET METRICS
# coverage metrics - bedtools
echo "calculating coverage metrics"
bedtools coverage \
    -abam ${prefix}_primer_trimmed.bam -b merged_targets_5col.bed -d \
    > ${prefix}.covd
awk '{sum+=$7}END{m=(sum/NR); b=m*0.2; print m, b}' ${prefix}.covd \
    > ${prefix}_covd.tmp 2> ${prefix}_coverage_metrics.log

# ontarget metrics - picard tools
# make intervals file for CollectTargetedPcrMetrics
samtools view -H ${prefix}_primer_trimmed.bam > ${prefix}_header.txt
cat ${prefix}_header.txt merged_targets_5col.bed > ${prefix}_fullintervals
cat ${prefix}_header.txt merged_targets_5col.bed > ${prefix}_noprimerintervals

echo "calculating on-target metrics"
$picard CollectTargetedPcrMetrics I=${prefix}_primer_trimmed.bam \
    O=${prefix}_targetPCRmetrics.txt AI=${prefix}_fullintervals \
    TI=${prefix}_noprimerintervals R=$ref \
    PER_TARGET_COVERAGE=${prefix}_perTargetCov.txt \
    2> ${prefix}_ontarget_pcrmetrics.log

# VARIANT CALLING
# GATK
echo "variant calling with GATK"
$gatk  -T HaplotypeCaller -R $ref -I ${prefix}_primer_trimmed.bam \
    -stand_call_conf 20 -stand_emit_conf 20 -mbq 20 -L merged_targets_5col.bed \
    -o ${prefix}_gatkHC.vcf 2> ${prefix}_gatkHC.log
# LoFreq
echo "variant calling with LoFreq"
lofreq call -q 20 -Q 20 -m 30 -C 50 --call-indels \
    -f $ref ${prefix}_primer_trimmed.bam -l merged_targets_5col.bed \
    -o ${prefix}_lofreq.vcf 4> ${prefix}_lofreq.log

echo "completed"
done
