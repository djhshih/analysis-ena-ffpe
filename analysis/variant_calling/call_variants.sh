#!/bin/bash

## Citation: 
## https://gatk.broadinstitute.org/hc/en-us/articles/360035531132--How-to-Call-somatic-mutations-using-GATK4-Mutect2
## https://gatk.broadinstitute.org/hc/en-us/community/posts/12450796994459-Asking-for-advice-on-Mutect2-calling-in-somatic-but-amplicon-data
## https://gatk.broadinstitute.org/hc/en-us/community/posts/4410094938395-Mutect2-on-targeted-single-cell-data-from-tumor

set -euo pipefail

if (( $# < 1 )); then
	echo "usage: ${0##*/} <bam-file>"
	exit 1
fi

root_outdir="../../data/vcf"

bam_path=$1
ref_path="../../data/ref/Homo_sapiens_assembly38.fasta"
germline_resource="../../data/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz"
panel_of_normals="../../data/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz"
small_exac_common="../../data/gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz"
index_bundle="../../data/gatk-test-data/mutect2/Homo_sapiens_assembly38.index_bundle"

if [ ! -f $bam_path ]; then
    echo "$bam_path does not exist"
    exit 1
fi

bam_basename=${bam_path##*/}
sample_id=${bam_basename%%.*}

outdir="$root_outdir/$sample_id"
mkdir -p $outdir

echo "Inputs:"
echo "Germline Resource: $germline_resource"
echo "Panel of Normals: $panel_of_normals"
echo "Reference Genome: $ref_path"
echo "BAM: $bam_path"
echo -e "# -----------------------------------\n"


echo -e "\nRunning GATK Mutect2...\n"

gatk Mutect2 \
    -R $ref_path \
    -I $bam_path \
    -germline-resource $germline_resource \
    -pon $panel_of_normals \
    --f1r2-tar-gz "${outdir}/${sample_id}_f1r2.tar.gz" \
    -O "${outdir}/${sample_id}_unfiltered.vcf" \
    --disable-read-filter NotDuplicateReadFilter \
    --downsampling-stride 50 \
    --linked-de-bruijn-graph \
    --max-reads-per-alignment-start 0


echo -e "\nRunning GATK LearnReadOrientationModel...\n"

gatk LearnReadOrientationModel \
    -I "${outdir}/${sample_id}_f1r2.tar.gz" \
    -O "${outdir}/${sample_id}_read-orientation-model.tar.gz"


echo -e "\nRunning GATK GetPileupSummaries...\n"

gatk GetPileupSummaries \
    -I $bam_path \
    -V $small_exac_common \
    -L $small_exac_common \
    -O "${outdir}/${sample_id}_pileupsummaries.table"


echo -e "\nRunning GATK CalculateContamination...\n"

gatk CalculateContamination \
    -I "${outdir}/${sample_id}_pileupsummaries.table" \
    -tumor-segmentation "${outdir}/${sample_id}_segments.table" \
    -O "${outdir}/${sample_id}_contamination.table"


echo -e "\nRunning GATK FilterMutectCalls...\n"

gatk FilterMutectCalls \
    -V "${outdir}/${sample_id}_unfiltered.vcf" \
    -R $ref_path \
    --tumor-segmentation "${outdir}/${sample_id}_segments.table" \
    --contamination-table "${outdir}/${sample_id}_contamination.table" \
    --ob-priors "${outdir}/${sample_id}_read-orientation-model.tar.gz" \
    -O "${outdir}/${sample_id}_filtered.vcf"


echo -e "\nRunning GATK FilterAlignmentArtifacts...\n"

gatk FilterAlignmentArtifacts \
    -R $ref_path \
    -V "${outdir}/${sample_id}_filtered.vcf" \
    -I $bam_path \
    --bwa-mem-index-image $index_bundle \
    -O "${outdir}/${sample_id}_alignment_artifacts_filtered.vcf"


echo -e "\n##########################################################"
echo -e "Finished Variant Calling for ${sample_id}\n"
