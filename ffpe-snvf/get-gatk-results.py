#!/usr/bin/env python
import polars as pl
from pathlib import Path
import sys
from tqdm import tqdm

repo_root = Path("..")

## Local Dependencies
sys.path.append(str(repo_root / "common-ffpe-snvf/python"))
from common import read_variants

## Functions
def process_dataset(dataset: str, variant_set: str, vcf_ext: str = "vcf.gz") -> None:
	"""Read the variant and filter columns from VCF and write 
	a table for evaluating the gatk orientation bias mixture model.

	Parameters:
		dataset: Dataset identifier (e.g. 'FFX').
		variant_set: Variant set folder name (e.g. 'mutect2-tn_filtered_pass-orientation').
	"""
	print(f"Processing... Dataset: {dataset} | Variant Set: {variant_set}")

	outdir_root = (repo_root / "ffpe-snvf" / dataset / variant_set / "gatk-obmm")
	vcf_paths = list((repo_root / "vcf" / dataset / variant_set).glob(f"*/*.{vcf_ext}"))

	tags = [f"Sample_B83_00{n}" for n in range(17, 32+1)] + ["FFPE"]
	vcf_paths = [
		path for path in vcf_paths 
		if any(tag in str(path) for tag in tags)
	]

	for path in tqdm(vcf_paths):
		sample_name = path.parent.name
		outdir = (outdir_root / sample_name)
		outdir.mkdir(exist_ok=True, parents=True)

		gatk_res = read_variants(path, columns=["#CHROM", "POS", "REF", "ALT", "FILTER"])
		gatk_res.write_csv((outdir / sample_name).with_suffix(".gatk-obmm.tsv"), separator="\t")

## Get GATK OBMM results
# process_dataset(dataset = "PRJEB8754", variant_set = "dup-unmarked_filtered_pass-orient-pos-sb")
process_dataset(dataset = "PRJEB44073", variant_set = "filtered_pass-orientation")
process_dataset(dataset = "SRP044740", variant_set = "filtered_pass-orientation")
process_dataset(dataset = "SRP065941", variant_set = "filtered_pass-orientation")
