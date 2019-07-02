## primerclip

#### primerclip v0.3.9beta primer trimming tool.

#### Swift Biosciences Inc. 2017

Primerclip™ is an alignment-based primer trimming tool designed to trim
primer sequences for Swift Biosciences Accel-Amplicon™ panels. The motivation
for designing an alignment-based primer trimming tool to increase speed.
Primerclip uses an algorithm based on genomic intervals of the aligned reads
rather than the sequence matching approach necessarily used by most trimming
tools to trim adapter sequences. Trimming based on alignment position allows
primerclip to run in significantly less time than sequence-based trim tools,
particularly as the size of the targeted panel increases.

If you have questions or would like additional support, please contact
Swift Technical Support at:

email: TechSupport@swiftbiosci.com
phone: 734 330 2568

#### Accel-Amplicon Trimming and Analysis Files for each Panel
Please visit [Swift Biosciences Accel-Amplicon™ trim and analysis files](https://swiftbiosci.com/protected-content/protected-content_amplicon-bed-files/)
to download the files needed for trimming primers with primerclip and performing
targeted variant calling.

__CHANGELOG__

190205 uncommented alternate chromosome parser to enable new arbitrary chromosome
parsing code.

190130 merged "extchrnames" branch to master to fix issue with some chromosome
names failing to parse. Primerclip should now parse arbitrary chromosome names
successfully.

__INSTALL PRE-COMPILED BINARY:__

This binary is compiled for linux on x86_64 (Ubuntu 16.04).

Path to the pre-compiled binary:

    .stack-work/install/x86_64-linux/lts-11.0/8.2.2/bin/primerclip

To install, copy the "primerclip" file into a folder on your PATH, e.g /usr/local/bin

    cp .stack-work/install/x86_64-linux/lts-11.0/8.2.2/bin/primerclip /usr/local/bin

Test that the binary is accessible by running

    which primerclip

If the path to the binary is returned, the binary should be ready to use.
If no path is returned, please double-check the location of the binary,
ensure that the folder with the binary is on your PATH:

    echo $PATH

and ensure that execute permissions are set for the binary:

    chmod a+x /path/to/primerclip

NOTE: although the pre-compiled binary is statically linked, it will still
      require that glibc version 2.19 is installed for the binary to run
      successfully. If you don't have this version of glibc
      ("ldd --version" to check), you will need to build from source
      (see instructions below).

__USAGE:__

To trim paired-end alignments:

    primerclip masterfile.txt alignmentfile.sam outputfilename.sam

To trim single-end alignments:

    primerclip -s masterfile.txt alignmentfile.sam outputfilename.sam

OPTIONAL USE OF BEDPE format for primer coordinates:

    primerclip -b primer_coords_bedpe.bed alignmentfile.sam outputfilename.sam

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

3. Run primerclip to generate a primer-trimmed SAM file.


Primerclip takes three arguments and one optional switch:

1. alignment input file in SAM format (see above if you're not sure how to get a
    SAM file from your FASTQ files) [alignmentfile.sam]

2. a text "master" input file containing the target and primer
    coordinates (as well as other information) for each amplicon in the panel
    used to create the library. [masterfile.sam]
    Master files for each Swift Accel-Amplicon panel are available at
    swiftbiosci.com. Primerclip extracts the primer genomic coordinates from
    the master file to use for trimming primers from the aligned reads.

   optionally, the "-b" or "--bedpe" switch can be included and a BEDPE primer
   coordinates file can be used as the second argument if your primer coordinates
   are not in master file format.


3. the desired name of the output file containing the primer-trimmed
    alignments in SAM format. [outputfilename.sam]

The primerclip output file is in SAM format; it is common to convert The
output SAM file to sorted and indexed BAM format using samtools or
picard-tools. The sorted, indexed BAM can then be used for downstream analysis
steps.

__BUILDING FROM SOURCE:__

This project is organized to use the haskell-stack build tool http://www.haskellstack.org.
To build the project from source, please follow the instructions on the
stack site for installing the stack build tool, then clone this repository
and run the following inside the primerclip project folder:

    stack build
    stack install

(stack install installs binary in ~/.local/bin, which can be added to your PATH)

The primerclip binary has been tested on the following operating systems:

* Ubuntu 14.04
* Ubuntu 16.04
* macOS Sierra

