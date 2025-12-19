#!/usr/bin/env python3

# Prepare json input files for WDL fastq_aligned_paired workflow 

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
ref_root = os.path.abspath('../../../ref')
ref_fname = 'Homo_sapiens_assembly38.fasta'

infname = '../../../../annot/PRJEB8754/sample-info_matched-ff-ffpe_on-pat-id-sample-type.tsv'
rg_root = os.path.abspath('../../rg')
fq_root = os.path.abspath('../../fq')
out_dir = 'inputs'


if not os.path.exists(out_dir):
    os.makedirs(out_dir)

fastq_pd = pd.read_csv(infname, sep='\t')

# ----

samples = []
for _, x in fastq_pd.groupby('run_accession'):
    sample_name = x['sample_name'].iloc[0]
    alias = x['sample_alias'].iloc[0]
    fastqs_r1 = return_path_if_exists(f"{fq_root}/{alias}/{sample_name}_1.fastq.gz")
    fastqs_r2 = return_path_if_exists(f"{fq_root}/{alias}/{sample_name}_2.fastq.gz")
    samples.append(
        {
            'sample_name': f"{sample_name}",
            # 'alias': f"{alias}",
            'fastqs_r1': fastqs_r1,
            'fastqs_r2': fastqs_r2,
        }
    )

ref_fpath = os.path.join(ref_root, ref_fname)

# base wdl input
base = {
    'fastq_align_paired.fastq_bwa_mem_paired.ref_fasta': return_path_if_exists(ref_fpath),
    'fastq_align_paired.fastq_bwa_mem_paired.ref_fasta_amb': return_path_if_exists(ref_fpath + '.64.amb'),
    'fastq_align_paired.fastq_bwa_mem_paired.ref_fasta_ann': return_path_if_exists(ref_fpath + '.64.ann'),
    'fastq_align_paired.fastq_bwa_mem_paired.ref_fasta_alt': return_path_if_exists(ref_fpath + '.64.alt'),
    'fastq_align_paired.fastq_bwa_mem_paired.ref_fasta_pac': return_path_if_exists(ref_fpath + '.64.pac'),
    'fastq_align_paired.fastq_bwa_mem_paired.ref_fasta_bwt': return_path_if_exists(ref_fpath + '.64.bwt'),
    'fastq_align_paired.fastq_bwa_mem_paired.ref_fasta_sa': return_path_if_exists(ref_fpath + '.64.sa'),
    'fastq_align_paired.fastq_bwa_mem_paired.cpu': 4,
    'fastq_align_paired.fastq_bwa_mem_paired.memory_gb': 12,
    'fastq_align_paired.fastq_bwa_mem_paired.preemptible_tries': 1,
    'fastq_align_paired.bam_sort_coord.cpu': 4,
    'fastq_align_paired.bam_sort_coord.memory_gb': 4,
    'fastq_align_paired.bam_sort_coord.preemptible_tries': 1,
}

# write wdl input json file for each sample
for x in samples:
    out = base.copy()
    out['fastq_align_paired.sample_name'] = x['sample_name']
    out['fastq_align_paired.fastq_bwa_mem_paired.fastqs_r1'] = [x['fastqs_r1']]
    out['fastq_align_paired.fastq_bwa_mem_paired.fastqs_r2'] = [x['fastqs_r2']]
    out['fastq_align_paired.fastq_bwa_mem_paired.rg_header'] = return_path_if_exists(f"{rg_root}/{x["sample_name"]}.rg.txt")
    out_path = os.path.join(out_dir, x['sample_name'] + '.inputs')
    with open(out_path, 'w') as outf:
        outf.write(json.dumps(out, indent=True, sort_keys=True))

