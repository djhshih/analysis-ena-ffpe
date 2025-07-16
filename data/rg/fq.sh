#!/bin/bash

# Extract the first entry from the fastq files

mkdir -p fq

for x in ../fq/*/*.fastq.gz; do
	fname=${x##*/}
	sid=${fname%%_*}
	zcat $x | head -4 > fq/${sid}.fq
done
