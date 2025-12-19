#!/usr/bin/env python3

# Prepare json input files for WDL mutect2 workflow 

import os
import json
import pandas as pd

def return_path_if_exists(path: str, abs=True) -> str:
	if os.path.exists(path):
		return os.path.abspath(path) if abs else path
	else:
		raise FileNotFoundError(f"File not found: {path}")

# ---

# absolute path is required in the json input files for WDL
ref_root = '../../../data/ref'
ref_fname = 'Homo_sapiens_assembly38'
vcf_root = "../../../data/gatk-best-practices/somatic-hg38"
bam_root = "../../../data/PRJEB8754/bam"
bundle_root = "../../../data/gatk-test-data/mutect2"
intervals_root = "../../../data/misc"
infname = '../../../annot/PRJEB8754/sample-info_matched-ff-ffpe_on-pat-id-sample-type.tsv'
gatk_path = "/home/moyukh/miniconda3/envs/ffpe-bench/share/gatk4-4.6.2.0-0/gatk-package-4.6.2.0-local.jar"

out_dir = 'inputs'

if not os.path.exists(out_dir):
	os.makedirs(out_dir)

ref_fpath = os.path.abspath(os.path.join(ref_root, ref_fname))
vcf_path = os.path.abspath(vcf_root)
bundle_path = os.path.abspath(bundle_root)
intervals_path = os.path.abspath(intervals_root)

pheno = pd.read_csv(infname, sep='\t')


# base wdl input
base = {
	'bam_variant_mutect2.run_funcotator': False,
	'bam_variant_mutect2.intervals': return_path_if_exists(os.path.join(ref_root, 'wgs_calling_regions.hg38.interval_list')),
	'bam_variant_mutect2.ref_fasta': return_path_if_exists(ref_fpath + '.fasta'),
	'bam_variant_mutect2.ref_fai': return_path_if_exists(ref_fpath + '.fasta.fai'),
	'bam_variant_mutect2.ref_dict': return_path_if_exists(ref_fpath + '.dict'),
	'bam_variant_mutect2.pon': return_path_if_exists(os.path.join(vcf_path, '1000g_pon.hg38.vcf.gz')),
	'bam_variant_mutect2.pon_idx': return_path_if_exists(os.path.join(vcf_path, '1000g_pon.hg38.vcf.gz.tbi')),
	'bam_variant_mutect2.gnomad': return_path_if_exists(os.path.join(vcf_path, 'af-only-gnomad.hg38.vcf.gz')),
	'bam_variant_mutect2.gnomad_idx': return_path_if_exists(os.path.join(vcf_path, 'af-only-gnomad.hg38.vcf.gz.tbi')),
	'bam_variant_mutect2.variants_for_contamination':  return_path_if_exists(os.path.join(vcf_path, 'small_exac_common_3.hg38.vcf.gz')),
	'bam_variant_mutect2.variants_for_contamination_idx': return_path_if_exists(os.path.join(vcf_path, 'small_exac_common_3.hg38.vcf.gz.tbi')),
	# 'bam_variant_mutect2.realignment_index_bundle': os.path.join(bundle_path, 'Homo_sapiens_assembly38.index_bundle'),
	'bam_variant_mutect2.m2_extra_args': '--disable-read-filter NotDuplicateReadFilter --downsampling-stride 50 --linked-de-bruijn-graph --max-reads-per-alignment-start 0', #  --max-reads-per-alignment-start 500 --dont-use-soft-clipped-bases --annotations-to-exclude StrandBiasBySample --annotations-to-exclude ReadPosRankSumTest
	'bam_variant_mutect2.scatter_count': 4,
	'bam_variant_mutect2.gatk_docker': 'broadinstitute/gatk:4.6.2.0',
	'bam_variant_mutect2.gatk_override': return_path_if_exists(gatk_path),
	'bam_variant_mutect2.bam_mutect2.mem': 4,
	'bam_variant_mutect2.run_orientation_bias_mixture_model_filter': True,
}

# Citations: https://gatk.broadinstitute.org/hc/en-us/community/posts/12450796994459-Asking-for-advice-on-Mutect2-calling-in-somatic-but-amplicon-data
# ----

for i in range(pheno.shape[0]):
	
	sample_name = pheno.loc[i, "sample_name"]
	
	bam = f"{bam_root}/{sample_name}.bam"
	bai = f"{bam_root}/{sample_name}.bai"
	
	if not os.path.exists(bam):
		raise FileNotFoundError(f"{bam} does not exist")
	if not os.path.exists(bai):
		raise FileNotFoundError(f"{bai} does not exist")

	# write wdl input json file for each sample
	out = base.copy()
	out['bam_variant_mutect2.tumor_bam'] = bam
	out['bam_variant_mutect2.tumor_bai'] = bai
	out_path = os.path.join(out_dir, f"{sample_name}.inputs")
	with open(out_path, 'w') as outf:
		outf.write(json.dumps(out, indent=True, sort_keys=True))

	print(f"Prepared inputs for: {sample_name}")

print("Done.")

