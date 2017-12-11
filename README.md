Swift Accel-Amplicon Primerclip v0.1 primer trimming tool.
Swift Biosciences Inc. 2017

If you have questions or would like to provide feedback, please contact us:

    TechSupport@swiftbiosci.com

Primerclip is an alignment-based primer trimming tool designed for trimming
primer sequences from aligned reads sequenced from a targeted high-throughput
sequencing library generated from PCR-amplified library molecules. The motivation
for designing an alignment-based primer trimming tool is speed: primerclip uses
an algorithm based on genomic intervals of the aligned reads rather than the
sequence matching approach necessarily used by most trimming tools to trim
adapter sequence. Trimming based on alignment position allows primerclip to run
in significantly less time than sequence-based trim tools, particularly as the
size of the targeted panel increases.

TODO: add table of run-times on control samples comparing cutadapt to primerclip

INSTALL PRE-COMPILED BINARY:

This binary is compiled for linux on x86_64 (Ubuntu 16.04).

To install, copy the "primerclip" file into a folder on your PATH, e.g /usr/local/bin

Test that the binary is accessible by running
    which primerclip

If the path to the binary is returned, the binary should be ready to use.
If no path is returned, please double-check the location of the binary,
ensure that the folder with the binary is on your path (echo $PATH),
and ensure that execute permissions are set for the binary:

  chmod a+x /path/to/primerclip

NOTE: although the pre-compiled binary is statically linked, it will still
      require that glibc version 2.19 is installed for the binary to run
      successfully. If you don't have this version of glibc
      ("ldd --version" to check), you will need to build from source
      (see instructions below).

USAGE:

    primerclip alignmentfile.sam masterfile.txt outputfilename.sam

This primer trimming tool is designed to be used after the Illumina adapter
trimming and alignment steps:

    1. Trim Illumina adapters (if present), trim low-quality bases, and
       remove short reads (input should be one or two FASTQ files)

    2. Align to reference genome with aligner of your choice (we recommend
       bwa mem). Input is trimmed FASTQ file(s), output should be in SAM
       format. NOTE: specifying an output filename with the ".sam" suffix
       should result in a SAM output file from bwa mem; if your input file is
       in BAM foramt, you can use the samtools command below to convert BAM to SAM:

            samtools view -h yourbamfile.bam > yoursamfile.sam

    3. You should now have an alignment file in SAM format which can be used as
       input to the primerclip tool.

 Primerclip takes three arguments:

   1. alignment file in SAM format (see above if you're not sure how to get a
      SAM file from your FASTQ files)

   2. a text "master" file containing the target and primer
      coordinates (as well as other information) for each amplicon in the panel
      used to create the library.
      Master files for each Swift Accel-Amplicon panel are available at
      swiftbiosci.com. Primerclip extracts the primer genomic coordinates from
      the master file to use for trimming primers from the aligned reads.

   3. the desired name of the output file containing the primer-trimmed
      alignments in SAM format.

The primerclip output file is in SAM format; it is usually a good idea to
convert the output SAM file to sorted and indexed BAM format using samtools or
picard-tools. The sorted, indexed BAM can then be used for downstream analysis
steps.

BUILDING FROM SOURCE:

This project is organized to use the haskell-stack build tool (www.haskellstack.org).
To build the source, follow the instructions on the stack site for installing
the stack build tool, then clone this repository and run the following:

1. stack init
2. stack build
3. stack install (installs the binary in ~/.local/bin)

The primerclip binary has been tested on the following operating systems:
    - Ubuntu 14.04
    - Ubuntu 16.04
    - macOS Sierra
