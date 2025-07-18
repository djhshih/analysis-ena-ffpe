#!/usr/bin/env python3

# Prepare json input files for WDL mutect2 workflow 

import os
import json
import pandas as pd
import numpy as np

# ---

# absolute path is required in the json input files for WDL
ref_root = '../../data/ref'
ref_fname = 'Homo_sapiens_assembly38'
vcf_root = "../../data/gatk-best-practices/somatic-hg38"
bam_root = "../../data/bam"
bundle_root = "../../data/gatk-test-data/mutect2"
intervals_root = "../../data/misc"
infname = '../../annot/sample-info_matched-ff-ffpe_on-pat-id.tsv'
gatk_path = "/home/moyukh/miniconda3/envs/ena-ffpe/share/gatk4-4.6.2.0-0/gatk-package-4.6.2.0-local.jar"

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
	'bam_variant_mutect2.intervals': os.path.join(intervals_path, 'standard_chromosomes.list'),
	'bam_variant_mutect2.ref_fasta': ref_fpath + '.fasta',
	'bam_variant_mutect2.ref_fai': ref_fpath + '.fasta.fai',
	'bam_variant_mutect2.ref_dict': ref_fpath + '.dict',
	'bam_variant_mutect2.pon': os.path.join(vcf_path, '1000g_pon.hg38.vcf.gz'),
	'bam_variant_mutect2.pon_idx': os.path.join(vcf_path, '1000g_pon.hg38.vcf.gz.tbi'),
	'bam_variant_mutect2.gnomad': os.path.join(vcf_path, 'af-only-gnomad.hg38.vcf.gz'),
	'bam_variant_mutect2.gnomad_idx': os.path.join(vcf_path, 'af-only-gnomad.hg38.vcf.gz.tbi'),
	'bam_variant_mutect2.variants_for_contamination':  os.path.join(vcf_path, 'small_exac_common_3.hg38.vcf.gz'),
	'bam_variant_mutect2.variants_for_contamination_idx': os.path.join(vcf_path, 'small_exac_common_3.hg38.vcf.gz.tbi'),
	'bam_variant_mutect2.realignment_index_bundle': os.path.join(bundle_path, 'Homo_sapiens_assembly38.index_bundle'),
	'bam_variant_mutect2.m2_extra_args': '--disable-read-filter NotDuplicateReadFilter --downsampling-stride 100 --linked-de-bruijn-graph --max-reads-per-alignment-start 0', # --dont-use-soft-clipped-bases --max-reads-per-alignment-start 500
	'bam_variant_mutect2.scatter_count': 4,
	'bam_variant_mutect2.gatk_docker': 'broadinstitute/gatk:4.6.2.0',
	'bam_variant_mutect2.gatk_override': gatk_path,
	'bam_variant_mutect2.bam_mutect2.mem': 8,
}

# Citations: https://gatk.broadinstitute.org/hc/en-us/community/posts/12450796994459-Asking-for-advice-on-Mutect2-calling-in-somatic-but-amplicon-data
# ----


for i in range(pheno.shape[0]):
	
	alias = pheno.loc[i, "sample_alias"]
	run_accession = pheno.loc[i, "run_accession"]
	
	bam = f"{bam_root}/{alias}_{run_accession}.bam"
	bai = f"{bam_root}/{alias}_{run_accession}.bai"
	
	if not os.path.exists(bam):
		raise FileNotFoundError(f"{bam} does not exist")
	if not os.path.exists(bai):
		raise FileNotFoundError(f"{bai} does not exist")

	# write wdl input json file for each sample
	out = base.copy()
	out['bam_variant_mutect2.tumor_bam'] = bam
	out['bam_variant_mutect2.tumor_bai'] = bai
	out_path = os.path.join(out_dir, f"{alias}_{run_accession}.inputs")
	with open(out_path, 'w') as outf:
		outf.write(json.dumps(out, indent=True, sort_keys=True))

	print(f"Prepared inputs for: {alias}")

print("Done.")

