# De Bruijn Graph-based Assembler

You can find the complete description of the practical work [here]( 
https://docs.google.com/document/d/1P4v3bHbSurD7RXA-ldVwmtNKGvWsnBae51RMGye_KLs/edit?usp=sharing).

## Introduction

The goal of this practical work is to assemble the genome of the Enterovirus A71. This genome is particularly interesting because it is very short: 7408 nucleotides, linear, and non-segmented. The fastq file you have was generated using the ART program [Huang 2011] with the following command: art_illumina -i eva71.fna -ef -l 100 -f 20 -o eva71 -ir 0 -dr 0 -ir2 0 -dr2 0 -na -qL 41 -rs 1539952693 The reads have maximum quality (41) and contain no insertions. Only the reads corresponding to the 5’ -> 3’ strands are provided here.

In the debruijn-tp/data/ folder, you will find:

    eva71.fna: genome of the virus of interest
    eva71_plus_perfect.fq: reads

## Installing Dependencies

You will use the Python libraries networkx, pytest, and pylint:

```
pip3 install --user networkx pytest pylint pytest-cov
```

## Usage

Run the Python3 program named debruijn.py in the debruijn/ folder. It will take as arguments: -i single-end fastq file -k kmer size (optional - default 21) -o output file with contigs.


