# Gap2Seq

Gap2Seq is a program for filling gaps in scaffolds produced by genome assembly
tools using short read data such as reads produced by Illumina sequencing.

Since version 3.0, it can also genotype insertions from variant calls produced
by variant calling tools.

## Reference

L. Salmela, K. Sahlin, V. Mäkinen, and A.I. Tomescu: Gap filling as
exact path length problem. In Proc. RECOMB 2015, LNBI 9029, Springer
2015, pp. 281-292.

L. Salmela, A.I. Tomescu: Safely filling gaps with partial solutions
common to all solutions. In Proc. WABI 2016, LNBI 9838, Springer
2016, xiv, short abstract.

R. Walve, L. Salmela, V. Mäkinen: Variant Genotyping with Gap Filling.
(Submitted)

## Requirements

Gap2Seq has been tested on systems running Linux on a X86_64 architecture.
Gap2Seq uses GATB library (http://gatb-core.gforge.inria.fr/index.html) for the
de Bruijn graph implementation and htslib (http://www.htslib.org) for reading
alignments for read filtering. The libraries are included in the Gap2Seq
package.

Compiling Gap2Seq requires gcc version 4.5 or newer and cmake.

## Installation

Unpack the Gap2Seq package.
For compiling Gap2Seq run

```
mkdir build;  cd build;  cmake ..;  make
```

The main script Gap2Seq and the binaries Gap2Seq-core, GapMerger, GapCutter and
ReadFilter can then be found in the build directory.

## Usage

```
Gap2Seq [parameters]

Required parameters:
-f, --filled <FASTA file>   output file for filled scaffolds
-r, --reads <FASTA/Q files> short reads, several files can be specified as a list separated by ','

-s, --scaffolds <FASTA/Q file>

-g, --gaps <FASTA/Q file>
-b, --bed <BED file>

Optional parameters:
-l, --libraries <lib conf>   

-t, --threads <int>         number of threads to use  [default 1]
-k <int>                    kmer length for DBG  [default 31]
--max-mem <float>           maximum memory usage of DP table computation in gigabytes (excluding DBG) [default 20]
--fuz <int>                 number of nucleotides to ignore on gap fringes  [default 10]
--dist-error <int>          maximum error in gap estimates  [default 500]
--solid <int>               threshold for solid k-mers for building the DBG [default 2]

-h, --help                  show this help message and exit
```

## Example

This example shows how to run Gap2Seq on the GAGE S. aureus data.

Download the GAGE data sets from

http://gage.cbcb.umd.edu/data/Staphylococcus_aureus/Data.original.tgz
http://gage.cbcb.umd.edu/data/Staphylococcus_aureus/Assembly.tgz

Unpack the data files.

### Regular

Run Gap2Seq (here we run it for the SGA scaffolds)

```
Gap2Seq --scaffolds Assembly/SGA/genome.scf.fasta --filled Assembly/SGA/genome.scf.fill.fasta --reads Data/original/frag_1.fastq,Data/original/frag_2.fastq,Data/original/shortjump_1.fastq,Data/original/shortjump_2.fastq
```

The filled scaffolds are then in the file Assembly/SGA/genome.scf.fill.fasta.

### With read filtering

Align, sort, and index the read libraries to the scaffolds with e.g. BWA MEM.

```
bwa index Assembly/SGA/genome.scf.fasta

bwa mem Assembly/SGA/genome.scf.fasta Data/original/frag_1.fastq Data/original/frag_2.fastq | samtools sort -O bam - > Data/original/frag.bam
samtools index Data/original/frag.bam

bwa mem Assembly/SGA/genome.scf.fasta Data/original/shortjump_1.fastq Data/original/shortjump_2.fastq | samtools sort -O bam - > Data/original/shortjump.bam
samtools index Data/original/shortjump.bam
```

Create a read library configuration file, a tab-separated list with a single
read library per line:
bam, read length, mean insert size, std dev, threshold

```
echo -e "Data/original/frag.bam\t101\t180\t18\t40" > libraries.txt
echo -e "Data/original/shortjump.bam\t37\t3500\t350\t40" >> libraries.txt
```

Run Gap2Seq.

```
Gap2Seq --scaffolds Assembly/SGA/genome.scf.fasta --filled Assembly/SGA/genome.scf.fill.fasta --libraries libraries.txt
```

## Changelog
### Version 3.0

Gap2Seq.sh is replaced with a Python script, which accepts gaps/scaffolds in
FASTA/FASTQ and VCF formats and reads in FASTA/FASTQ and SAM/BAM formats.

Optional per-gap read filtering when run with new script.

Gap2Seq binary is renamed Gap2Seq-core and can still be used instead of the new
script.

Flanks of length between k and k+fuz are now used rather than the hard limit of
k+fuz.

Switched to GATB-core 1.2.2.

### Version 2.0

Parallelization is now on gap level when run with the Gap2Seq.sh script.

Safe bases inserted into gaps are outputted in upper case and unsafe
bases are outputted in lower case.

### Version 1.0

Optimized version of the algorithm.

Output now indicates on which gaps search was aborted because of the
memory limit.

### Version 0.3

Reorganized parallel execution.

### Version 0.2

Proper synchronization for access to the memuse hash table.

Switched to GATB 1.0.5.

The maximum memory limitation option is now total for all threads.
This is then divided evenly to all threads.

Memory usage tracking now includes all major data structures excluding
the DBG.
