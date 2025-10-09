# %%
#!/usr/bin/env python3
import polars as pl
import os

# Read in sample info sheet
metadata = pl.read_csv("sample-info_unprocessed.tsv", separator="\t")


# %%
sample_annot = (
	metadata
	.with_columns(
		(pl.col("sample_alias").str.split("_").alias("sample_alias_split")),
    )
	.with_columns(
        pl.col("sample_alias_split").list.get(0).str.split("-").list.get(0).alias("case_id"),
        pl.col("sample_alias_split").list.get(1).alias("sample_type"),
        pl.col("sample_alias_split").list.get(2).alias("preservation"),
		(pl.col("sample_alias") + pl.lit("_") + pl.col("run_accession")).alias("sample_name")
    )
	.drop("sample_alias_split")
	.select(['sample_name','sample_alias','case_id','sample_type','preservation', 'run_accession','sample_accession','experiment_accession','study_accession','instrument_model','library_name','fastq_ftp',])
	.sort("case_id", "preservation")
)

# Save to file
sample_annot.write_csv("sample-info_stage2.tsv", separator="\t")
print("Saved table with all sample information to: sample-info_stage1.tsv\n")

# %%
# Obtain sample counts
# 10 matched FFPE FF samples as mentioned in the pubication
print("Counting FFPE and Frozen samples available")
matched_ffpe_ff = (
    sample_annot
	.group_by(["case_id", "sample_type"], maintain_order=True)
	.agg([
		(pl.col("preservation") == "FFPE").sum().alias("FFPE_count"),
		(pl.col("preservation") == "Frozen").sum().alias("FF_count"),
	])
)

matched_ffpe_ff.write_csv("sample_preservation_count.tsv", separator="\t")
print("\tCount table saved to: sample_preservation_count.tsv\n")

# %%
## Obtain counts were there is atleast 1 FFPE and 1 FF sample
matched_ffpe_ff = matched_ffpe_ff.filter((pl.col("FFPE_count") > 0) & (pl.col("FF_count") > 0))

# Filter just based on the Patients who have 1 FFPE and 1 FF sample. 
# This also keeps the primary FFPE samples along with Metastasis samples
sample_info_matched_ffpe_ff_pat = (
	sample_annot
 	.filter(pl.col("case_id").is_in(matched_ffpe_ff.get_column("case_id").implode()))
)

# Filtering for matched FFPE and FF withing the same patient
# Keeping both Primary and Metastasis
sample_info_matched_ffpe_ff_pat.write_csv("sample-info_matched-ff-ffpe_on-pat-id.tsv", separator="\t")
print("Filtered Sample info with the criteria of at least 1 FFPE and 1 FF sample from each patient: sample-info_matched-ff-ffpe_on-pat-id.tsv")


# %%

## Filtering based on Patients who have 1 FFPE and 1 FF sample within the same sample type. 
## This only keeps the 10 matched FFPE-FF metastasis samples as mentioned in the publications
sample_info_matched_ffpe_ff_pat_stype = sample_annot.join(matched_ffpe_ff, on=["case_id", "sample_type"], how="semi")
sample_info_matched_ffpe_ff_pat_stype.write_csv("sample-info_matched-ff-ffpe_on-pat-id-sample-type.tsv", separator="\t")
print("Filtered Sample info with the criteria of atleast 1 FFPE and 1 FF sample from each patient within same sample type: sample-info_matched-ff-ffpe_on-pat-id-sample-type.tsv")

# %%
## Create a bash script to dowload the fastqs
fastq_links = (
	sample_annot
 	.with_columns(pl.col("fastq_ftp").str.split(";"))
  	.explode("fastq_ftp")
)

fastq_links.select("fastq_ftp").write_csv("fq_ftp_links.txt", include_header=False, separator="\t")
print("Saved fastq ftp links to: fq_ftp_links.txt")

wget = ["#!/bin/bash"]

for i in range(fastq_links.shape[0]):
    link = f"ftp://{fastq_links[i, "fastq_ftp"]}"
    sample_name = f"{fastq_links[i, "sample_name"]}"
    
    wget.append(f'mkdir -p {sample_name} && wget "ftp://{link}" -O {sample_name}/{sample_name}.fastq.gz')


fastq_get_path = "../../data/PRJEB8754/fq/fastq_ftp_download.sh"
with open(fastq_get_path, "w") as file:
    for line in wget:
        file.write(f"{line}\n")

print(f"Created a bash script to obtain fastqs from the EBI ftp at : {fastq_get_path}")



