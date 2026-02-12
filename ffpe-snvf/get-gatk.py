#!/usr/bin/env python
import polars as pl
import os
import glob
from tqdm import tqdm

def read_variants(path:str, columns: list = ["#CHROM", "POS", "REF", "ALT", "FILTER"]) -> pl.DataFrame:
	variants = (
		pl.read_csv(path, separator="\t", comment_prefix="##", infer_schema_length=1000, columns=columns)
		.rename(lambda x: x.lstrip("#").lower())
		.with_columns(pl.col("alt").str.split(","))
		.explode("alt")
	)
	return variants

## VCFs from each dataset to use for getting gatk obmm results
dset_vcf = {
	"PRJEB8754" : "filtered_pass-orient-pos-sb-ad-blacklist-macni_dup-unmarked",
	"PRJEB44073" : "filtered_pass-orientation-dp20-blacklist-macni",
	"SRP044740" : "filtered_pass-orientation-dp20-blacklist-macni",
	"SRP065941" : "filtered_pass-orientation-dp20-blacklist"
}

## Process each datasets
for i, dset in enumerate(dset_vcf.keys()):
	
	var_set = dset_vcf[dset]
	print(f"Processing Dataset: {dset} | {var_set}")

	vcf_paths = glob.glob(f"../vcf/{dset}/{var_set}/*/*.vcf")

	for path in tqdm(vcf_paths):
		sample_name = path.split("/")[-2]

		if "froz" in sample_name.lower():
			continue

		variants = read_variants(path)

		outdir = f"{dset}/{var_set.removeprefix("vcf_")}/gatk-obmm/{sample_name}"
		os.makedirs(outdir, exist_ok=True)

		variants.write_csv(f"{outdir}/{sample_name}.gatk-obmm.tsv", separator="\t")

