#!/usr/bin/env bash
set -euox pipefail

N_SAMPLES=40
AVG_THRESH=15
DATASET="SRP044740"
CHR_SIZES="../Homo_sapiens_assembly38.chrom.sizes"

zcat "${DATASET}.SiteDepth.gz" | \
	awk -v n="$N_SAMPLES" -v thresh="$AVG_THRESH" '$3 >= thresh*n {print $1"\t"$2"\t"$2+1}' | \
	bedtools merge -i - > "${DATASET}_pileup-avg-dp15.bed"

bedtools slop -b 200 -g "$CHR_SIZES" -i "${DATASET}_pileup-avg-dp15.bed" | \
	bedtools merge -i - > "${DATASET}_pileup-avg-dp15-padded.bed"
