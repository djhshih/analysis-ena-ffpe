#!/bin/bash

# Extract the first entry from the fastq files

mkdir -p fq_head

for x in ../fq/*/*.fastq.gz; do
	fname=${x##*/}
	sid=${fname%%_*}
	zcat $x | head -4 > fq_head/${sid}.fq
done
