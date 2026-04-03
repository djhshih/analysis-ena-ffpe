#!/usr/bin/env bash

set -euo pipefail

# Define the function
filter_vcf() {
    # Check if the required arguments are provided
    if [[ $# -lt 3 ]]; then
        echo "Usage: filter_vcf <indir_root> <outdir_root> <dataset_name>"
        return 1
    fi

    local indir_root="$1"
    local outdir_root="$2"
    local dataset_name="$3"
    
    # Define target directory based on the output directory name
    local target_indir="../../germline-filter/$dataset_name/$(basename "$outdir_root")"

    echo "Input directory: $indir_root"
    echo "Output directory: $outdir_root"
    echo "Target directory: $(realpath $target_indir)"

    local i=1
    for vcf in "$indir_root"/*/*.vcf.gz; do
        
        # Safety check: skip if no files match the glob
        [[ -e "$vcf" ]] || continue

        echo -e "\t$i. Filtering: $vcf"

        local sample_name=$(basename "$(dirname "$vcf")")
        local target_list="$target_indir/$sample_name/$sample_name.tsv"

        if [[ -f "$target_list" ]]; then

            local outdir="$outdir_root/$sample_name"
            mkdir -p "$outdir"
            local outpath="$outdir/$sample_name.vcf.gz"

            bcftools view -T "$target_list" "$vcf" -O z -o "$outpath"
            bcftools index -t "$outpath"

        else
            echo -e "\t   [WARNING] $target_list does not exist. Skipping..."
        fi

        ((i++))

    done

    echo -e "Outputs saved to: $outdir_root\n"
}

# --- Execution ---

filter_vcf "filtered_pass-orientation-exome-blacklist" "filtered_pass-orientation-exome-blacklist-macni" "PRJEB44073"


