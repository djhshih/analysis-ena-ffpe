import os
import glob
import polars as pl

## Functions
def check_res_existence(path: str, abs = False, halt: bool = False) -> str | None:
	if os.path.exists(path):
		return os.path.abspath(path) if abs else path
	else:
		if halt:
			raise FileNotFoundError(f"{path} does not exist.")
		else:
			print(f"'{path}' does not exist. skipping...")
			return None


def filter_micr_artifacts(ffpe_filter: str, dataset: str, variant_sets: dict, sample: str, micr_artifacts: pl.DataFrame, res_suffix: str = ".snv") -> pl.DataFrame:

	print(f"\tFiltering MICR artifacts from {ffpe_filter.upper()} results...")
	ffpe_filter_res_path = check_res_existence(f"{dataset}/{variant_sets[dataset]}/{ffpe_filter.lower()}/{sample}/{sample}.{ffpe_filter.lower()}{res_suffix}")
 
	if ffpe_filter_res_path is not None:
		ffpe_filter_res = pl.read_csv(ffpe_filter_res_path, separator="\t")
		ffpe_filter_res_micr_filtered = ffpe_filter_res.join(micr_artifacts, on =["chrom", "pos", "ref", "alt"], how="anti")

		output_dir = f"{dataset}/{variant_sets[dataset].replace("filtered", "filtered_micr")}/{ffpe_filter.lower()}/{sample}"
		os.makedirs(output_dir, exist_ok=True)
		output_path = f"{output_dir}/{sample}.{ffpe_filter.lower()}{res_suffix}"

		ffpe_filter_res_micr_filtered.write_csv(output_path, separator="\t")
		print(f"\t\tMICR filtered {ffpe_filter.upper()} results written to: {output_path}")

		variants_before = ffpe_filter_res.height
		artifact_variants = ffpe_filter_res.height - ffpe_filter_res_micr_filtered.height
		variants_after = ffpe_filter_res_micr_filtered.height
		print(f"\t\t{artifact_variants} MICR artifacts removed from {variants_before} variants ({(artifact_variants) / variants_before * 100:.2f}%).")
		print(f"\t\t{variants_after} variants remain after MICR filtering.")

		summary = pl.DataFrame({
			"dataset": dataset,
			"sample": sample,
			"ffpe_filter" : f"{ffpe_filter.lower()}",
			"total_variants": variants_before,
			"micr_artifacts_removed": artifact_variants,
			"variants_after_micr_filtering": variants_after,
			"percent_removed": artifact_variants / variants_before * 100,
		})

		return summary

	elif ffpe_filter_res_path is None:
		return pl.DataFrame()


## Setup
variant_sets = {
	"PRJEB8754" : "vcf_pass-orient-pos-sb_ad_filtered",
	"PRJEB44073" : "vcf_filtered_pass_orientation",
	"SRP044740" : "vcf_filtered_pass_orientation",
	"SRP065941" : "vcf_filtered_pass_orientation",
}

datasets = list(variant_sets.keys())

## Filter MICR artifacts from each dataset

for dataset in datasets:

	microsec_paths = glob.glob(f"{dataset}/{variant_sets[dataset]}/microsec/outputs/*/*.microsec.tsv")

	micr_filtering_summaries = []

	for i, path in enumerate(microsec_paths):
		sample = os.path.basename(path).replace(".microsec.tsv", "")
		print(f"\n{i+1}. Processing... Dataset: {dataset}, Sample: {sample}")

		microsec = pl.read_csv(path, separator="\t", infer_schema_length=1000).rename({"Chr": "chrom"}).rename(lambda x: x.lower())
		microsec_artifacts = microsec.filter(pl.col("msec_filter_all") == "Artifact suspicious")

		mobsnvf_summary = filter_micr_artifacts("mobsnvf", dataset, variant_sets, sample, microsec_artifacts)
		micr_filtering_summaries.append(mobsnvf_summary)
	
		sobdetector_summary = filter_micr_artifacts("sobdetector", dataset, variant_sets, sample, microsec_artifacts)
		micr_filtering_summaries.append(sobdetector_summary)

		vafsnvf_summary = filter_micr_artifacts("vafsnvf", dataset, variant_sets, sample, microsec_artifacts)
		micr_filtering_summaries.append(vafsnvf_summary)
	
		ideafix_summary = filter_micr_artifacts("ideafix", dataset, variant_sets, sample, microsec_artifacts, res_suffix=".tsv")
		micr_filtering_summaries.append(ideafix_summary)

	micr_filtering_summaries_df = pl.concat(micr_filtering_summaries, how="diagonal_relaxed")
	micr_filtering_summaries_df.write_csv(f"{dataset}/{variant_sets[dataset].replace("filtered", "filtered_micr")}/micr_filtering_summary.tsv", separator="\t")
	


