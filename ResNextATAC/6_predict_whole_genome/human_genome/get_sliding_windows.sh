#!/usr/bin/env bash

# take human chr1 as an example
GENOME=/path/to/Genome/hg38/GRCh38.p14.genome.fa
## make sliding windows: 500bp window size, 100bp step size
bedtools makewindows -b Human_chr1.bed -w 500 -s 100 > Region.chr1.slide.Human.bed
bedtools getfasta -tab -fi $GENOME -bed Region.chr1.slide.Human.bed -fo Region.chr1.slide.Human.input.txt
