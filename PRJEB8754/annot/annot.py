#!/usr/bin/env python3
import polars as pl
import os
import requests
import xml.etree.ElementTree as ET

# Read in sample info sheet
sample_info = pl.read_csv("sample-info_stage0.tsv", separator="\t")
sample_info

## The current sample info does not have sufficient information about the samples.
## They are obtained from the EBI ENA API
annotation = {
	"sample_accession" : [],
	"sample_title" : [],
	"sample_alias" : [],
	"center_name" : []
}

print("Obtaining sample information from EBI ENA API for:")
for sample_accession in sample_info.get_column("sample_accession"):
    print(f"\t{sample_accession}")
    
    response = requests.get(f"https://www.ebi.ac.uk/ena/browser/api/xml/{sample_accession}?includeLinks=true")
    response.raise_for_status()
    root = ET.fromstring(response.content)
    
    annotation["sample_accession"].append(sample_accession)
    annotation["sample_title"].append(root.find("SAMPLE").find("TITLE").text)
    annotation["sample_alias"].append(root.find("SAMPLE").get("alias"))
    annotation["center_name"].append(root.find("SAMPLE").get("center_name"))
    
extra_info = pl.DataFrame(annotation)

# Join information obtained with current sample info table
sample_info_extra = sample_info.join(extra_info, on = "sample_accession")

print("Done\n")

# Refactor columns
sample_info_new = (
    sample_info_extra
    .with_columns(
        pl.col("sample_alias").str.split("_").alias("sample_alias_split")
    )
    .with_columns(
        pl.col("sample_alias_split").list.get(0).alias("inferred_id"),
        pl.col("sample_alias_split").list.get(1).alias("sample_type"),
        pl.col("sample_alias_split").list.get(2).alias("preservation")
    )
    .drop("sample_alias_split") # Optionally drop the intermediate list column
    .select(['sample_title',  'inferred_id', 'sample_type', 'preservation', 'run_accession', 'sample_accession', 'experiment_accession', 'study_accession',  'sample_alias', 'center_name', 'tax_id', 'scientific_name', 'fastq_ftp', 'sra_ftp', 'bam_ftp'])
    .sort(["inferred_id", "preservation"])
)

sample_info_new.write_csv("sample-info_stage1.tsv", separator="\t")
# sample_info_new

# Ssave to file
sample_info_new = pl.read_csv("sample-info_stage1.tsv", separator="\t")
print("Saved table with additional sample information to: sample-info_stage1.tsv\n")



# sample_info_new.filter(pl.col("sample_type") == "Prim")
# sample_info_new.filter(pl.col("inferred_id") == "Pat05")

# Obtain sample counts
# 10 matched FFPE FF samples as mentioned in the pubication
print("Counting FFPE and Frozen samples available")
matched_ffpe_ff = (
    sample_info_new
	.group_by(["inferred_id", "sample_type"], maintain_order=True)
	.agg([
		(pl.col("preservation") == "FFPE").sum().alias("FFPE_count"),
		(pl.col("preservation") == "Frozen").sum().alias("FF_count"),
	])
)

matched_ffpe_ff.write_csv("sample_preservation_count.tsv", separator="\t")
print("\tCount table saved to: sample_preservation_count.tsv\n")

## Obtain counts were there is atleast 1 FFPE and 1 FF sample
matched_ffpe_ff = matched_ffpe_ff.filter((pl.col("FFPE_count") > 0) & (pl.col("FF_count") > 0))
matched_ffpe_ff

# Filter just based on the Patients who have 1 FFPE and 1 FF sample. This also keeps the primary FFPE samples along with Metastasis samples
sample_info_matched_ffpe_ff_pat = (
	sample_info_new
 	.filter(
      pl.col("inferred_id").is_in(matched_ffpe_ff.get_column("inferred_id").implode())
    )
)

# Keeping both Primary and Metastasis
sample_info_matched_ffpe_ff_pat.write_csv("sample-info_matched-ff-ffpe_on-pat-id.tsv", separator="\t")
print("Filtered Sample info with the criteria of atleast 1 FFPE and 1 FF sample from each patient: sample-info_matched-ff-ffpe_on-pat-id.tsv")
# sample_info_matched_ffpe_ff_pat

## Filtering based on Patients who have 1 FFPE and 1 FF sample within the same sample type. 
## This only keeps the 10 matched FFPE-FF metastasis samples as mentioned in the publications
sample_info_matched_ffpe_ff_pat_stype = sample_info_new.join(matched_ffpe_ff, on=["inferred_id", "sample_type"], how="semi")
sample_info_matched_ffpe_ff_pat_stype.write_csv("sample-info_matched-ff-ffpe_on-pat-id-sample-type.tsv", separator="\t")
print("Filtered Sample info with the criteria of atleast 1 FFPE and 1 FF sample from each patient within same sample type: sample-info_matched-ff-ffpe_on-pat-id-sample-type.tsv")

# sample_info_matched_ffpe_ff_pat_stype

## Create a bash script to dowload the fastqs
fastq_links = (
	sample_info
 	.with_columns(pl.col("fastq_ftp").str.split(";"))
	.select("fastq_ftp")
  	.explode("fastq_ftp")
	.with_columns([
		(pl.lit("ftp://") + pl.col("fastq_ftp")).alias("fastq_ftp")
	])
)

fastq_links.write_csv("fq_ftp_links.txt", include_header=False, separator="\t")
print("Saved fastq ftp links to: fq_ftp_links.txt")
# fastq_links

wget = ["#!/bin/bash"]

for link in fastq_links.get_column("fastq_ftp"):
    basename = os.path.basename(link)
    dir_name = f"{link.split("/")[-2]}"
    
    wget.append(f"mkdir -p {dir_name} && wget {link} -O {dir_name}/{basename}")


fastq_get_path = "../data/fq/fastq_ftp_download.sh"
with open(fastq_get_path, "w") as file:
    for line in wget:
        file.write(f"{line}\n")

print(f"Created a bash script to obtain fastqs from the EBI ftp at : {fastq_get_path}")


