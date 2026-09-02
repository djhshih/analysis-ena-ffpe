# %%
import os
os.environ["POLARS_MAX_THREADS"] = "1"

import polars as pl
import subprocess
from pathlib import Path
from scipy.stats import false_discovery_control
import numpy as np



# %% [markdown]
# ## Setup

# %%
project_root = Path("..").resolve()
vcf_root = project_root / "vcf"
oxog_snvf_root = project_root / "oxog-snvf"

# %% [markdown]
# ## Functions

# %%
def select_oxog_context(df: pl.DataFrame, invert: bool = False) -> pl.DataFrame:
	condition = ((pl.col("ref") == "C") & (pl.col("alt") == "A")) | (
		(pl.col("ref") == "G") & (pl.col("alt") == "T")
	)

	if invert:
		return df.filter(~condition)
	else:
		return df.filter(condition)


def annotate_adaptive_fdr(df: pl.DataFrame) -> pl.DataFrame:

	df = df.with_columns(pl.col("FOBP").cast(pl.Float64))
	
	# Change scores of 0 to machine epsilon
	df = df.with_columns(
		pl.when(pl.col("FOBP") == 0)
		.then(np.finfo(float).eps)
		.otherwise(pl.col("FOBP"))
		.alias("score")
	)

	# Adjust MOBSNVF FOBP using FDR BH method
	df = df.with_columns(
		pl.col("score")
		.map_batches(
			lambda s: false_discovery_control(s.to_numpy(), method="bh"),
			return_dtype=pl.Float64,
		)
		.alias("q")
	)

	# Multiply Q values with their ordinal rank
	df = df.with_columns(
		(pl.col("q").rank(method="ordinal") * pl.col("q")).alias("adaptive_fdr")
	)

	df = df.drop("score")

	return df

# %%
def process(dataset, source_variant_set, oxog_variant_set, threshold = 1e-13):

	vcf_dir = vcf_root / dataset /  source_variant_set
	oxog_snvf_dir = oxog_snvf_root / dataset / oxog_variant_set / "mobsnvf"

	print(f"Dataset: {dataset} | Source Variant Set: {source_variant_set} | OXOG Variant Set: {oxog_variant_set}")

	vcf_paths = list(vcf_dir.glob("*/*.vcf.gz"))

	for i, path in enumerate(vcf_paths, start=1):
		sample_name = path.name.split(".")[0]

		oxog_res_path = oxog_snvf_dir / sample_name / f"{sample_name}.mobsnvf.snv"

		oxog_res = pl.read_csv(oxog_res_path, separator="\t", infer_schema_length=1000).pipe(select_oxog_context)

		if oxog_variant_set != source_variant_set:
			vcf = pl.read_csv(path, separator="\t", comment_prefix="##", columns=["#CHROM", "POS", "REF", "ALT"]).rename(lambda x: x.lstrip("#").lower())
			oxog_res = oxog_res.join(vcf, how="semi", on=["chrom", "pos", "ref", "alt"])

		oxog_artifacts = oxog_res.pipe(annotate_adaptive_fdr).filter(pl.col("adaptive_fdr") >= threshold)

		vcf_outdir = vcf_root / dataset / f"{source_variant_set}-oxogsnvf" / sample_name
		vcf_outdir.mkdir(parents=True, exist_ok=True)

		out_vcf_path = vcf_outdir / f"{sample_name}.vcf.gz"

		artifacts_path = vcf_outdir / f"{sample_name}.oxog_artifacts.tsv"
		oxog_artifacts.write_csv(artifacts_path, separator="\t", include_header=False)

		print(f"{i}. Processing sample: {sample_name}")
		if not oxog_artifacts.is_empty():
			subprocess.run(
				[
					"bcftools", "view",
					"-T", f"^{artifacts_path}",
					"-Oz", "-o", str(out_vcf_path),
					str(path),
				],
				check=True,
			)
		else:
			subprocess.run(
				[
					"bcftools", "view",
					"-Oz", "-o", str(out_vcf_path),
					str(path),
				],
				check=True,
			)

		subprocess.run(["bcftools", "index", "-t", str(out_vcf_path)], check=True)


# %% [markdown]
# ## Run Filters

# %%
process(
	dataset = "SRP065941", 
	source_variant_set = "filtered_pass-orientation-exome-blacklist", 
	oxog_variant_set = "filtered_pass-orientation-exome", 
	threshold = 1e-4
)

# %%
process(
	dataset = "PRJEB44073", 
	source_variant_set = "filtered_pass-orientation-exome-blacklist-macni", 
	oxog_variant_set = "filtered_pass-orientation-exome", 
	threshold = 1e-4
)

# %%
process(
	dataset = "SRP044740", 
	source_variant_set = "filtered_pass-alignment-exome-blacklist-macni", 
	oxog_variant_set = "filtered_pass-orientation-exome-blacklist", 
	threshold = 1e-12
)


