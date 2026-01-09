#!/usr/bin/env python
import polars as pl
import glob
import os


def get_ffpe_snvf_paths(dataset: str, variant_set: str) -> list:
	"""
	Returns the path for each FFPE SNVF model's
	results for a specified variant set and dataset
	"""

	paths = sorted(
		glob.glob(f"{dataset}/{variant_set}/*/*/*.mobsnvf.snv") +
		glob.glob(f"{dataset}/{variant_set}/*/*/*.vafsnvf.snv") +
		glob.glob(f"{dataset}/{variant_set}/*/*/*.sobdetector.snv") +
		glob.glob(f"{dataset}/{variant_set}/*/*/*.microsec.tsv") +
		glob.glob(f"{dataset}/{variant_set}/*/*/*.ideafix.tsv") +
		glob.glob(f"{dataset}/{variant_set}/*/*/*.gatk-obmm.tsv") +
		glob.glob(f"{dataset}/{variant_set}/*/*/*.ffpolish.tsv")
	)

	return paths


## These are strictly not necessary as the source and target variant sets for all the datasets share the same structure
## However this can be useful if this naming convention ever breaks
## Links dataset to variant source sets to filter from
dset_vset = {
	"PRJEB44073" : "vcf_filtered_pass_orientation",
	"SRP044740" : "vcf_filtered_pass_orientation",
	"SRP065941" : "vcf_filtered_pass_orientation"
}

# ## Target DP filtered vcf used to subset variants from DP unfiltered variant set
dset_target_vcf = {
	"PRJEB44073" : "vcf_filtered_pass-orientation-dp10",
	"SRP044740" : "vcf_filtered_pass-orientation-dp10",
	"SRP065941" : "vcf_filtered_pass-orientation-dp10"
}

## Process each datasets
for i, dset in enumerate(dset_vset.keys()):
	
	var_set = dset_vset[dset]
	print(f"Processing Dataset: {dset} | {var_set}")

	## Get paths to all the SNVF results
	paths = get_ffpe_snvf_paths(dset, var_set)

	## Subset each sample
	for i, path in enumerate(paths):

		model = os.path.basename(path).split(".")[-2]
		sample_name = os.path.basename(path).split(".")[0]
		new_varset = dset_target_vcf[dset]

		print(f"\t{i+1}. Processing sample: {sample_name}")

		## SNV set from model output to filter
		snvf = pl.read_csv(path, separator="\t", infer_schema_length=1000)

		## Read in the filtered VCF, to subset the snvf results with
		vcf_path = (f"../vcf/{dset}/{new_varset}/{sample_name}/{sample_name}.vcf")

		vcf_df = (
			pl.read_csv(
				vcf_path, 
				separator="\t", 
				comment_prefix="##", 
				null_values=".", 
				columns=["#CHROM", "POS", "REF", "ALT"], 
				infer_schema_length=1000
			)
			.rename(lambda x: x.lstrip("#").lower())
			.with_columns(pl.col("alt").str.split(",")).explode("alt")
		)

		if model == "microsec":
			subset = snvf.join(vcf_df, left_on=["Chr", "Pos", "Ref", "Alt"], right_on=["chrom", "pos", "ref", "alt"], how="semi")
		else:
			subset = snvf.join(vcf_df, on=["chrom", "pos", "ref", "alt"], how="semi")

		
		# Write output
		new_snvf_name = new_varset.removeprefix("vcf").lstrip("_").lstrip("-")
		output_dir = f"{dset}/{new_snvf_name}/{model}/{sample_name}"
		os.makedirs(output_dir, exist_ok=True)

		if model in ["mobsnvf", "vafsnvf", "sobdetector"]:
			output_path = f"{output_dir}/{sample_name}.{model}.snv"
		else:
			output_path = f"{output_dir}/{sample_name}.{model}.tsv"
		subset.write_csv(output_path, separator="\t")
	
		print(f"\tWritten subsetted data to {output_path}\n")
			

print("All files processed.")



