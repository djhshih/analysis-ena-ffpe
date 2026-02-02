import polars as pl
import glob
import os
import vcformer
from tqdm import tqdm
import seaborn as sns

datasets = ["PRJEB44073", "SRP044740"]

dset_vset = {
    "SRP044740": "vcf_filtered_pass-orientation-dp10-blacklist",
    "PRJEB44073":  "vcf_filtered_pass-orientation-dp10-blacklist"
}

for dataset in datasets:
	print(f"Processing dataset: {dataset}")
	paths = sorted(glob.glob(f"../vcf/{dataset}/{dset_vset[dataset]}/*/*.vcf"))

	all_variants = []

	for path in tqdm(paths):

		sample_name = path.split("/")[-2]
		
		vcf = (
			vcformer.read_vcf_as_polars(path,  sample_fields=["AD", "AF", "DP"], info_fields=["POPAF"])
			.rename(lambda x : x.replace(f"{sample_name}.", "").lower())
			.explode("alts", "af", "filters", "popaf")
			.with_columns(
				pl.col("ad").list.get(0).alias("ref_ad"),
				pl.col("ad").list.get(1).alias("alt_ad"),
				pl.lit(sample_name).alias("sample_name"),
				(10 ** -pl.col("popaf")).alias("popaf")
			)
			.drop("id", "qual", "phased", "ad")
			.rename({"alts": "alt"})
		)

		all_variants.append(vcf)

	outdir = dataset
	os.makedirs(outdir, exist_ok=True)

	all_variants_df = pl.concat(all_variants)
	all_variants_df.write_parquet(f"{dataset}/{dataset}_all-variants_{dset_vset[dataset].removeprefix("vcf_")}.parquet")


# %% [markdown]
# ### Explore

ax = sns.histplot(vcf["popaf"], bins=100)
ax.set_yscale("log")

vcf.filter(vcf["popaf"] > 0.02/100)


