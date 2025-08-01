#!/usr/bin/python3
import polars as pl
import pysam
import pandas as pd
import numpy as np
import os
import glob

outdir_root="../../ffpe-snvf/pass-n-orientation_ad_filtered/microsec/inputs"

## Create Mutation Info for MicroSEC

def import_formatted_vcf(vcf_path) -> pl.DataFrame:
    """
    Reads a VCF file and returns a Polars DataFrame formatted for MicroSEC analysis.

    The function:
    - Reads a VCF file as a tab-delimited table, ignoring header lines starting with '##'.
    - Drops all contigs, only retaining Chr 1-22, X, and Y
    - Selects the relevant columns (#CHROM, POS, REF, ALT).
    - Adds placeholder columns required by MicroSEC: Sample, Mut_type, SimpleRepeat_TRF, and Neighborhood_sequence.
    - Renames columns to standardized names: Chr, Pos, Ref, Alt.
    - Splits multi-allelic ALT entries into separate rows, resulting in a biallelic representation for each variant.

    Parameters
    ----------
    vcf_path : str
        Path to the VCF file to be imported.

    Returns
    -------
    pl.DataFrame
        A Polars DataFrame with columns: Sample, Mut_type, Chr, Pos, Ref, Alt, SimpleRepeat_TRF, Neighborhood_sequence.
        Each row represents a single biallelic variant.
    """
    allowed_chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    vcf = (
		pl.read_csv(vcf_path, separator="\t", comment_prefix="##")
		.filter(pl.col("#CHROM").is_in(allowed_chroms))
		.select(["#CHROM", "POS", "REF", "ALT"])
		.with_columns([
			pl.lit("-").alias("Sample"),
			pl.lit("-").alias("Mut_type"),
			pl.lit("-").alias("SimpleRepeat_TRF"),
			pl.lit("-").alias("Neighborhood_sequence")
		])
		.rename({
			"#CHROM" : "Chr",
			"POS" : "Pos", 
			"REF": "Ref", 
			"ALT" : "Alt"
		})
		.select(["Sample", "Mut_type", "Chr", "Pos", "Ref", "Alt", "SimpleRepeat_TRF", "Neighborhood_sequence"])
  		.with_columns(
			pl.col("Alt").str.split(",")
		)
    	.explode("Alt")
	)
    return vcf


lookup_table = pl.read_csv("../../annot/sample-info_matched-ff-ffpe-only.tsv", separator="\t")

vcf_paths = sorted(glob.glob("../../data/vcf_pass-n-orientation_ad_filtered/*/*.vcf"))
bam_paths = sorted(glob.glob(f"../../data/bam/*.bam"))

lookup_table = lookup_table.with_columns([
	pl.col("run_accession").map_elements(
		lambda id: [path for path in vcf_paths if id in path][0],
		return_dtype=pl.String
	).alias("vcf_path"),
  
 	pl.col("run_accession").map_elements(
		lambda name : [path for path in bam_paths if name in path][0],
		return_dtype=pl.String
	).alias("bam_path")
])

genome = pysam.FastaFile("../../data/ref/Homo_sapiens_assembly38.fasta")

lookup_table_ffpe = lookup_table.filter(pl.col("preservation") == "FFPE")

for i in range(lookup_table_ffpe.shape[0]):
    alias = lookup_table_ffpe[i, "sample_alias"]
    run_acc = lookup_table_ffpe[i, "run_accession"]
    sample_name = f"{alias}_{run_acc}"
    
    mut_info = import_formatted_vcf(lookup_table_ffpe[i, "vcf_path"]).with_columns(pl.lit(sample_name).alias("Sample"))
    
    outdir = f"{outdir_root}/mut_info"
    os.makedirs(outdir, exist_ok=True)
    
    for j in range(mut_info.shape[0]):
        
        chr_length = genome.get_reference_length(mut_info[j, "Chr"])
  
        left_flank = genome.fetch(
		mut_info[j, "Chr"],
		mut_info[j, "Pos"] - (21 if mut_info[j, "Pos"] > 20 else mut_info[j, "Pos"]),
		mut_info[j, "Pos"] - 1
		)
        
        right_flank = genome.fetch(
		mut_info[j, "Chr"],
		mut_info[j, "Pos"],
		mut_info[j, "Pos"] + (20 if (chr_length - mut_info[j, "Pos"] > 20) else (chr_length - mut_info[j, "Pos"]))
		)
        
        mut_info[j, "Neighborhood_sequence"] = (left_flank + mut_info[j, "Alt"] + right_flank).upper()
    
    snv_mask = mut_info["Ref"].str.len_chars() == mut_info["Alt"].str.len_chars()
    deletion_mask = (mut_info["Ref"].str.len_chars() > 1) & (mut_info["Alt"].str.len_chars() == 1)
    insertion_mask = (mut_info["Ref"].str.len_chars() == 1) & (mut_info["Alt"].str.len_chars() > 1)
    
    mut_info = (
		mut_info.with_columns([
			pl.when(snv_mask)
			.then(mut_info["Ref"].str.len_chars().cast(pl.Utf8) + "-snv")
			.when(deletion_mask)
			.then((mut_info["Ref"].str.len_chars() - 1).cast(pl.Utf8) + "-del")
			.when(insertion_mask)
			.then((mut_info["Alt"].str.len_chars() - 1).cast(pl.Utf8) + "-ins")
			.otherwise(mut_info["Mut_type"])
			.alias("Mut_type")
		])
	)
    
    # MicroSEC requires if variant to be at least 200 nucleotides away from the terminal ends of the chromosome
    mut_info = mut_info.with_columns(
		pl.col("Chr").map_elements(lambda chr: genome.get_reference_length(chr), return_dtype=pl.Int64).alias("chr_length")
	).filter(
		(pl.col("Pos") > 200) & (pl.col("Pos") < (pl.col("chr_length") - 200))
	).drop("chr_length")

    mut_info.write_csv(f"{outdir}/{sample_name}.microsec.mut-info.tsv", separator="\t")
    print(f"Finished creating Mut-Info for - {i + 1}. {sample_name}")


### Create Sample Info for MicroSEC

mut_info_paths = sorted(glob.glob(f"{outdir_root}/mut_info/*.microsec.mut-info.tsv"))
ref_path = '../../data/ref/Homo_sapiens_assembly38.fasta'
mut_info_suffix = ".microsec.mut-info.tsv"

sample_info = pd.DataFrame()

for i in range(len(mut_info_paths)):
    
    sample_name = os.path.basename(mut_info_paths[i]).split(mut_info_suffix)[0]
    run_acc = sample_name.split("_")[-1]
    
    lookup = lookup_table.filter(pl.col("preservation") == "FFPE").filter(pl.col("run_accession").str.contains(run_acc))
    
    sample_info.loc[i, "sample_name"] = sample_name
    sample_info.loc[i, "mutation information tsv file"] = mut_info_paths[i]
    sample_info.loc[i, "BAM file"] = lookup[0, "bam_path"]  #if i < len(bam_paths) else np.nan
    
    sample_info.loc[i, "read length"] = str(150) # majoriy read length in bam files
    sample_info.loc[i, "adapter sequence read 1"] = np.nan # DK if adapters are trimmed before alignment to bam
    sample_info.loc[i, "optional: adapter sequence read 2"] = np.nan # DK if reads were paired end
    sample_info.loc[i, "sample type: Human or Mouse"] = "hg38"
    sample_info.loc[i, "panel name"] = "TOP"
    sample_info.loc[i, "optional: reference genome fasta file"] = ref_path
    #sample_info.loc[i, "optional: simple repeat region bed file"] = np.nan

sample_info.to_csv(f"{outdir_root}/microsec.sample_info.tsv", sep="\t", index=False, header=False)
