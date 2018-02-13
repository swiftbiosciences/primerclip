#### Tool Versions:
1. bedtools: v2.17.0
2. samtools: 1.3.1
3. gatk: 3.6-0-g89b7209
4. lofreq: build-date-May 19 2015
5. picard: 1.129

#### Reference Genome:
human: Homo_sapiens_assembly19broad.fasta

#### File Descriptions:
1. nonmerged_targets.bed = bed file of targets designed for panel
2. merged_targets.bed = bedfile of overlapping targets  merged into single target
3. merged_targets_5col.bed = modified merged_targets bedfile for bedtool compatibility

4. *_adapter_trimmed.fq.gz = adapter trimmed fastq.gz
6. *_non_primer_trimmed.bam = aligned with bwa bam
8. *_primer_trimmed.bam = primer trimmed bam

9. *.covd = per base coverage metrics (per target) from bedtools
10. *_targetPCRmetrics.txt = pcr-related metrics for sample
11. *_perTargetCov.txt = coverage per target metrics for sample
12. *_gatkHC.vcf = gatk variant caller vcf
13. *_lofreq.vcf = lofreq variant caller vcf

#### Convert BAM to SAM
```bash
samtools view -h -o output.sam input.bam
```
