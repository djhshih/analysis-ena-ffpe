import polars as pl
import os
import glob


annot = pl.read_csv("sample_annotation_stage0.tsv", separator = "\t")

sample_annotation = (
	annot
	.select(["sample_alias", "run_alias", "experiment_alias", "sample_accession", "run_accession", "experiment_accession", "read_count", "scientific_name", "instrument_model", "fastq_ftp"])
	.sort("sample_alias")
	.with_columns(pl.col("sample_alias").str.split("_").alias("sample_alias_split"))
	.with_columns(
     	pl.col("sample_alias_split").list.get(-1).alias("preservation"),
		pl.col("sample_alias_split").list.get(0).str.split("").list.get(1).alias("case_id")
    )
	.with_columns(pl.when(pl.col("sample_alias").str.starts_with("N")).then(pl.lit("Normal")).otherwise(pl.lit("Tumor")).alias("sample_type"))
	.select(["sample_alias", "case_id", "sample_type", "preservation", "run_alias", "experiment_alias", "sample_accession", "run_accession", "experiment_accession", "read_count", "scientific_name", "instrument_model", "fastq_ftp"])
	.rename({"sample_alias":"sample_name"})
)


sample_annotation.write_csv("sample_annotation_stage1.tsv", separator = "\t")

fastq_annotation = (
	sample_annotation
	.select(["sample_name", "case_id", "sample_type", "preservation", "run_accession", "instrument_model", "fastq_ftp"])
	.with_columns(pl.col("fastq_ftp").str.split(";"))
	.explode("fastq_ftp")
)

fastq_annotation.write_csv("fastq_annotation.tsv", separator = "\t")

wget = ["#!/bin/bash"]

for i, link in enumerate(fastq_annotation.get_column("fastq_ftp")):
    basename = os.path.basename(link)
    run_accession = basename.split(".")[0].split("_")[0]
    sample_name = fastq_annotation.filter(pl.col("run_accession") == run_accession)[0, "sample_name"]
    
    dir_name = f"{link.split("/")[-2]}"
    
    wget.append(f"mkdir -p {sample_name} && wget -c ftp://{link} -O {sample_name}/{sample_name}.fastq.gz")


fq_outdir = "../../data/SRP065941/fq"
os.makedirs(fq_outdir, exist_ok=True)

fastq_get_path = f"{fq_outdir}/fastq_ftp_download.sh"
with open(fastq_get_path, "w") as file:
    for line in wget:
        file.write(f"{line}\n")

print(f"Created a bash script to obtain fastqs from the EBI ftp at : {fastq_get_path}")


