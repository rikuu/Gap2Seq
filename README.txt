Gap2Seq
Contact: leena.salmela@cs.helsinki.fi

--------
Overview
--------

Gap2Seq is a program for filling gaps in scaffolds produced by
genome assembly tools using short read data such as reads produced by
Illumina sequencing.

---------
Reference
---------

L. Salmela, K. Sahlin, V. Mäkinen, and A.I. Tomescu: Gap filling as
exact path length problem. In Proc. RECOMB 2015, LNBI 9029, Springer
2015, pp. 281-292.

L. Salmela, A.I. Tomescu: Safely filling gaps with partial solutions common to all solutions. (Submitted)

-------------------
System Requirements
-------------------

Gap2Seq has been tested on systems running Linux on a X86_64
architecture. Gap2Seq uses GATB library
(http://gatb-core.gforge.inria.fr/index.html) for the de Bruijn graph
implementation. The library is included in the Gap2Seq package.
Compiling Gap2Seq requires gcc version 4.5 or newer and cmake.

------------
Installation
------------

Unpack the Gap2Seq package.
For compiling Gap2Seq run

    mkdir build;  cd build;  cmake ..;  make

The main script Gap2Seq.sh and the binaries Gap2Seq, GapMerger and GapCutter can then be found in the build directory.

-----
Usage
-----

Gap2Seq.sh [parameters]

Required parameters:
-scaffolds <FASTA/Q file>    scaffolds to be gap filled
-filled <FASTA file>         output file for filled scaffolds
-reads <FASTA/Q files>       short reads, several files can be specified as a list separated by ','

Optional parameters:
-max-mem <float>             maximum memory usage of DP table computation in gigabytes (excluding DBG) [default 20]
-fuz <int>                   number of nucleotides to ignore on gap fringes  [default 10]
-dist-error <int>            maximum error in gap estimates  [default 500]
-solid                       threshold for solid k-mers for building the DBG [default 2]
-k                           kmer length for DBG  [default 31]
-all-upper		     If specified, all filled bases are in upper case.
-unique			     If specified, only gaps with a unique path of best length are filled.
-nb-cores                    number of cores to use [default 0 (all cores)]
-verbose                     verbosity level (currently does not affect much?)  [default 1]
-help                        display help about possible options

-------
Example
-------

This example shows how to run Gap2Seq on the GAGE S. aureus data.

Download the GAGE data sets from

http://gage.cbcb.umd.edu/data/Staphylococcus_aureus/Data.original.tgz
http://gage.cbcb.umd.edu/data/Staphylococcus_aureus/Assembly.tgz

Unpack the data files.

Run Gap2Seq (here we run it for the SGA scaffolds)

    Gap2Seq.sh -scaffolds Assembly/SGA/genome.scf.fasta -filled Assembly/SGA/genome.scf.fill.fasta -reads Data/original/frag_1.fastq,Data/original/frag_2.fastq,Data/original/shortjump_1.fastq,Data/original/shortjump_2.fastq

The filled scaffolds are then in the file Assembly/SGA/genome.scf.fill.fasta.

------------------
New in Version 2.0
------------------

Parallelization is now on gap level when run with the Gap2Seq.sh script.

Safe bases inserted into gaps are outputted in upper case and unsafe
bases are outputted in lower case.

------------------
New in Version 1.0
------------------

Optimized version of the algorithm.

Output now indicates on which gaps search was aborted because of the
memory limit.

------------------
New in Version 0.3
------------------

Reorganized parallel execution.

------------------
New in Version 0.2
------------------

Proper synchronization for access to the memuse hash table.

Switched to GATB 1.0.5.

The maximum memory limitation option is now total for all threads.
This is then divided evenly to all threads.

Memory usage tracking now includes all major data structures excluding
the DBG.
