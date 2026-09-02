import polars as pl
from pathlib import Path

project_root = Path("..").resolve()

annots = {
    "ENA Chong21": str(project_root / "annot/PRJEB44073/sample-info_stage3.tsv"),
    "ENA SRP044740": str(project_root / "annot/SRP044740/sample-info_stage2.tsv"),
    "ENA Oh15": str(project_root / "annot/SRP065941/sample_annotation_stage2_tumor-only.tsv")
}

accessions = {
    "ENA Chong21": "PRJEB44073",
    "ENA SRP044740": "SRP044740",
    "ENA Oh15": "SRP065941"
}

datasets = list(accessions.keys())

vcf_root = project_root / "vcf"
ff_variant_set = "filtered_pass-orientation-exome-blacklist-macni"

## Functions
def pool_ff(dataset:str, ff_variant_set:str):

	dset_acc = accessions[dataset]
	vcf_dir = vcf_root / dset_acc / ff_variant_set

	annot = pl.read_csv(annots[dataset], separator="\t")
	annot_ffpe = annot.filter(pl.col("preservation").str.to_lowercase().str.contains("ffpe"))
	annot_ff = annot.filter(pl.col("preservation").str.to_lowercase().str.contains("frozen"))
	ffpe_sample_names = annot_ffpe["sample_name"].unique()

	# return ffpe_sample_names

	outdir = Path(dset_acc) / "pooled_ff" / ff_variant_set
	outdir.mkdir(exist_ok=True, parents=True)

	for i, ffpe_sample_name in enumerate(ffpe_sample_names):

		case_id = annot_ffpe.filter(pl.col("sample_name") == ffpe_sample_name)[0, "case_id"]
		case_annot: pl.DataFrame = annot_ff.filter(pl.col("case_id") == case_id)
		n_ff: int = case_annot.height
		print(f"{i+1}. Pooling {n_ff} matched FF Variants from sample: {ffpe_sample_name}")

		pooled_ff: list = []
		for ff_sample_name in case_annot["sample_name"]:

			# return ff_sample_name
			vcf_path: Path = vcf_dir / ff_sample_name / f"{ff_sample_name}.vcf.gz"

			vcf_ff: pl.DataFrame = pl.read_csv(vcf_path, separator="\t", comment_prefix="##", columns=["#CHROM", "POS", "REF", "ALT"]).rename(lambda x : x.lstrip("#").lower())

			pooled_ff.append(vcf_ff)

		## Drop duplicate variants
		pooled_ff: pl.DataFrame = pl.concat(pooled_ff).unique()
		pooled_ff.write_csv(outdir / f"{ffpe_sample_name}.pooled-ff.tsv", separator="\t")


## Pool Fresh Frozens
pool_ff(dataset="ENA Chong21", ff_variant_set="filtered_pass-orientation-exome-blacklist-macni-oxogsnvf")

pool_ff(dataset="ENA Oh15", ff_variant_set="filtered_pass-orientation-exome-blacklist-oxogsnvf")

pool_ff(dataset="ENA SRP044740", ff_variant_set="filtered_pass-alignment-exome-blacklist-macni-oxogsnvf")


