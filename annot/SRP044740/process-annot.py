import polars as pl
import os
import glob

## SRP044740
# 
# Dataset with exomes of 13 FFPE breast tumor samples and 13 corresponding frozen samples.

sample_info = pl.read_csv("sample-info_stage0.tsv", separator="\t")
sample_info

new_sample_info = (
	sample_info
	.with_columns([
		pl.col("sample_alias").str.replace("BGI-FROZ", "").str.replace("BGI-FFPE", "").cast(int).alias("sample_number"),
		pl.col("sample_alias").str.replace("BGI-", "").str.replace(r"\d+$", "").alias("sample_type")
	])
	.sort("sample_number")
	.select(['study_title','sample_number','sample_type','run_accession','experiment_accession','sample_accession','scientific_name','sample_alias','instrument_model','fastq_md5','fastq_ftp'])
)

new_sample_info.write_csv("sample-info_stage1.tsv", separator="\t")

sample_count = (
	new_sample_info
	.group_by(["sample_number"])
	.agg([
		(pl.col("sample_type") == "FFPE").sum().alias("n_ffpe"),
		(pl.col("sample_type") == "FROZ").sum().alias("n_frozen")
	])
)

sample_count.write_csv("sample-count.tsv", separator="\t")

fastq_links = (
	new_sample_info
 	.with_columns(pl.col("fastq_ftp").str.split(";"))
	.select(["sample_alias", "fastq_ftp"])
  	.explode("fastq_ftp")
	.with_columns([
		(pl.lit("ftp://") + pl.col("fastq_ftp")).alias("fastq_ftp")
	])
)

fastq_links

fastq_links.select("fastq_ftp").write_csv("fq_ftp_links.txt", include_header=False, separator="\t")
print("Saved fastq ftp links to: fq_ftp_links.txt")


# fastq_links

wget = ["#!/bin/bash"]

for i in range(fastq_links.shape[0]):
    
    link = fastq_links[i, "fastq_ftp"]
    sample_name = fastq_links[i, "sample_alias"]
    
    basename = os.path.basename(link)
    read_number = basename.split("_")[1].split(".")[0]
    dir_name = f"{link.split("/")[-2]}"
    
    wget.append(f"mkdir -p {sample_name} && wget {link} -O {sample_name}/{sample_name}_{basename}")

fastq_dir = "../../data/SRP044740/fq"
os.makedirs(fastq_dir, exist_ok=True)

fastq_get_path = "../data/fq/fastq_ftp_download.sh"
with open(fastq_get_path, "w") as file:
    for line in wget:
        file.write(f"{line}\n")

print(f"Created a bash script to obtain fastqs from the EBI ftp at : {fastq_get_path}")

# parallel


script_dir = "../../data/SRP044740/fq/parallel_downloads"
os.makedirs(script_dir, exist_ok=True)

for i in range(fastq_links.shape[0]):
    
    link = fastq_links[i, "fastq_ftp"]
    sample_name = fastq_links[i, "sample_alias"]
    basename = os.path.basename(link)
    read_number = basename.split("_")[1].split(".")[0]
    dir_name = f"{link.split("/")[-2]}"
    
    wget = ["#!/bin/bash"]
    wget.append(f"mkdir -p {sample_name} && wget --tries=50 --retry-connrefused -O {sample_name}/{sample_name}_{basename} {link}")
    
    script_path = f"{script_dir}/{sample_name}_{basename}.sh"
    with open(script_path, "w") as file: 
        for line in wget:
        	file.write(f"{line}\n")
        



