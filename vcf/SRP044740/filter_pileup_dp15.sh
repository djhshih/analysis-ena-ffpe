#!/bin/env/bash

set -euo pipefail

outdir_root="filtered_pass-orientation-dp15pileup"
indir="filtered_pass-orientation/*/*.vcf.gz"
bed_dir="../../pandepth/SRP044740/individual_beds"

for vcf in $indir; do
    
    base=${vcf##*/}
    sample_name=${base%%.*}

    outdir=${outdir_root}/${sample_name}
    mkdir -p $outdir

    bed_path=${bed_dir}/${sample_name}.bed

    # ensure the bed file exists before attempting to filter
    if [[ ! -e "$bed_path" ]]; then
        echo "Warning: bed file not found: $bed_path" >&2
        echo "Skipping: $base"
        continue
    fi

    echo -e "\nFiltering: $base\n"

    bcftools view -T $bed_path -Oz -o ${outdir}/${sample_name}.vcf.gz $vcf
    bcftools index -t ${outdir}/${sample_name}.vcf.gz 

done

echo Done.
