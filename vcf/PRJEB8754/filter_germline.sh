#!/usr/bin/env bash

set -euo pipefail

target_indir="../../germline-filter/PRJEB8754/filtered_pass-orient-pos-sb-ad-blacklist-macni_dup-unmarked"
outdir_root=filtered_pass-orient-pos-sb-ad-blacklist-macni_dup-unmarked

for vcf in filtered_pass-orient-pos-sb-ad-blacklist_dup-unmarked/*/*.vcf; do

    echo $vcf

    sample_name=$(basename $(dirname $vcf))
    target_list=$target_indir/$sample_name/$sample_name.tsv

    if [[ -f $target_list ]]; then

        outdir=$outdir_root/$sample_name
        mkdir -p $outdir

        bcftools view -T $target_list $vcf -o $outdir/$sample_name.vcf

    else
        echo $target_list does not exist. 
        echo Find out why

    fi

done