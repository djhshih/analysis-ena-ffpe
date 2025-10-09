# %%
#!/usr/bin/python3
import polars as pl
import pysam
import pandas as pd
import numpy as np
import os
import glob


# %%

## Create Mutation Info for MicroSEC

def import_formatted_vcf(vcf_path: str, ct_only: bool = False) -> pl.DataFrame:
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
 
	if(ct_only):
		vcf = vcf.filter(
			((pl.col("Ref") == "C") & (pl.col("Alt") == "T")) |
			((pl.col("Ref") == "G") & (pl.col("Alt") == "A"))
		)

	return vcf


# %%
def make_mut_info(vcf_path: str, sample_name: str, genome: pysam.FastaFile, ct_only: bool = False) -> pl.DataFrame:
    
	mut_info = import_formatted_vcf(vcf_path, ct_only).with_columns(pl.lit(sample_name).alias("Sample"))
	
	# Define a function to be applied to each row-like struct
	# Used to retrieve neighbourhood sequence via map_elements which is vectorized for speed
	def get_neighborhood_sequence(row_struct: dict) -> str:
		chrom = row_struct["Chr"]
		pos = row_struct["Pos"]
		alt = row_struct["Alt"]

		try:
			chr_length = genome.get_reference_length(chrom)

			# Calculate start and end for fetching flanks, ensuring they are within bounds
			left_start = max(0, pos - 21)
			left_end = max(0, pos - 1)

			right_start = pos
			right_end = min(chr_length, pos + 20)

			left_flank = genome.fetch(chrom, left_start, left_end)
			right_flank = genome.fetch(chrom, right_start, right_end)

			return (left_flank + alt + right_flank).upper()

		except ValueError:
			# Pysam can raise ValueError for contigs not in the FASTA file
			return None

	# Create a struct of the columns needed for the neighborhood sequence calculation
	# and then apply the function.
	mut_info = mut_info.with_columns(
		pl.struct(["Chr", "Pos", "Alt"])
		.map_elements(get_neighborhood_sequence, return_dtype=pl.Utf8)
		.alias("Neighborhood_sequence")
	)
	
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
 
	return mut_info

# %%
## PRJEB8754

print("Making inputs for ENA PRJEB8754")

outdir_root = "../PRJEB8754/vcf_pass-orient-pos-sb_ad_filtered/microsec/inputs"
genome = pysam.FastaFile("../../data/ref/Homo_sapiens_assembly38.fasta")
lookup_table = pl.read_csv("../../annot/PRJEB8754/sample-info_matched-ff-ffpe-only.tsv", separator="\t")


vcf_paths = sorted(glob.glob("../../vcf/PRJEB8754/vcf_pass-orient-pos-sb_ad_filtered/*/*.vcf"))
bam_paths = sorted(glob.glob(f"../../data/PRJEB8754/bam/*.bam"))


lookup_table = lookup_table.with_columns([
	pl.col("run_accession").map_elements(
		lambda id: [os.path.abspath(path) for path in vcf_paths if id in path][0],
		return_dtype=pl.String
	).alias("vcf_path"),
  
 	pl.col("run_accession").map_elements(
		lambda name : [os.path.abspath(path) for path in bam_paths if name in path][0],
		return_dtype=pl.String
	).alias("bam_path")
])

lookup_table_ffpe = lookup_table.filter(pl.col("preservation") == "FFPE")

# %%
for i in range(lookup_table_ffpe.shape[0]):
	alias = lookup_table_ffpe[i, "sample_alias"]
	run_acc = lookup_table_ffpe[i, "run_accession"]
	sample_name = f"{alias}_{run_acc}"
	
 
	outdir = f"{outdir_root}/mut_info"
	os.makedirs(outdir, exist_ok=True)

	mut_info = make_mut_info(lookup_table_ffpe[i, "vcf_path"], sample_name, genome, ct_only = True)

	mut_info.write_csv(f"{outdir}/{sample_name}.microsec.mut-info.tsv", separator="\t")
	print(f"Finished creating Mut-Info for - {i + 1}. {sample_name}")


# %%
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
	sample_info.loc[i, "mutation information tsv file"] = os.path.abspath(mut_info_paths[i])
	sample_info.loc[i, "BAM file"] = lookup[0, "bam_path"]
	
	sample_info.loc[i, "read length"] = str(150) # majoriy read length in bam files
	sample_info.loc[i, "adapter sequence read 1"] = np.nan
	sample_info.loc[i, "optional: adapter sequence read 2"] = np.nan
	sample_info.loc[i, "sample type: Human or Mouse"] = "hg38"
	sample_info.loc[i, "panel name"] = "TOP"
	sample_info.loc[i, "optional: reference genome fasta file"] = ref_path
	#sample_info.loc[i, "optional: simple repeat region bed file"] = np.nan

sample_info.to_csv(f"{outdir_root}/microsec.sample_info.tsv", sep="\t", index=False, header=False)


# %%
## SRP044740

print("Making inputs for ENA SRP044740")

outdir_root = "../SRP044740/vcf_filtered_pass_orientation/microsec/inputs"
genome = pysam.FastaFile("../../data/ref/Homo_sapiens_assembly38.fasta")
lookup_table = pl.read_csv("../../annot/SRP044740/sample-info_stage1.tsv", separator="\t")


vcf_paths = sorted(glob.glob("../../vcf/SRP044740/vcf_filtered_pass_orientation/*/*.vcf"))
bam_paths = sorted(glob.glob(f"../../data/SRP044740/bam/*.bam"))


lookup_table = lookup_table.with_columns([
	pl.col("run_accession").map_elements(
		lambda id: [os.path.abspath(path) for path in vcf_paths if id in path][0],
		return_dtype=pl.String
	).alias("vcf_path"),
  
 	pl.col("run_accession").map_elements(
		lambda name : [os.path.abspath(path) for path in bam_paths if name in path][0],
		return_dtype=pl.String
	).alias("bam_path")
])

lookup_table_ffpe = lookup_table.filter(pl.col("sample_type") == "FFPE")
# lookup_table_ffpe

# %%
for i in range(lookup_table_ffpe.shape[0]):
	alias = lookup_table_ffpe[i, "sample_alias"]
	run_acc = lookup_table_ffpe[i, "run_accession"]
	sample_name = f"{alias}_{run_acc}"
	

	outdir = f"{outdir_root}/mut_info"
	os.makedirs(outdir, exist_ok=True)

	mut_info = make_mut_info(lookup_table_ffpe[i, "vcf_path"], sample_name, genome, ct_only = True)

	mut_info.write_csv(f"{outdir}/{sample_name}.microsec.mut-info.tsv", separator="\t")
	print(f"Finished creating Mut-Info for - {i + 1}. {sample_name}")

# %%
### Create Sample Info for MicroSEC

mut_info_paths = sorted(glob.glob(f"{outdir_root}/mut_info/*.microsec.mut-info.tsv"))
ref_path = os.path.abspath('../../data/ref/Homo_sapiens_assembly38.fasta')
mut_info_suffix = ".microsec.mut-info.tsv"

sample_info = pd.DataFrame()

for i in range(len(mut_info_paths)):
	
	sample_name = os.path.basename(mut_info_paths[i]).split(mut_info_suffix)[0]
	run_acc = sample_name.split("_")[-1]
	
	lookup = lookup_table.filter(pl.col("sample_type") == "FFPE").filter(pl.col("run_accession").str.contains(run_acc))

	sample_info.loc[i, "sample_name"] = sample_name
	sample_info.loc[i, "mutation information tsv file"] = os.path.abspath(mut_info_paths[i])
	sample_info.loc[i, "BAM file"] = lookup[0, "bam_path"]
	
	sample_info.loc[i, "read length"] = str(90) # majoriy read length in bam files
	sample_info.loc[i, "adapter sequence read 1"] = np.nan
	sample_info.loc[i, "optional: adapter sequence read 2"] = np.nan
	sample_info.loc[i, "sample type: Human or Mouse"] = "hg38"
	sample_info.loc[i, "panel name"] = "TOP"
	sample_info.loc[i, "optional: reference genome fasta file"] = ref_path
	#sample_info.loc[i, "optional: simple repeat region bed file"] = np.nan

sample_info.to_csv(f"{outdir_root}/microsec.sample_info.tsv", sep="\t", index=False, header=False)


# %%
## SRP065941

print("Making inputs for ENA SRP065941")

outdir_root = "../SRP065941/vcf_filtered_pass_orientation/microsec/inputs"
genome = pysam.FastaFile("../../data/ref/Homo_sapiens_assembly38.fasta")
lookup_table = pl.read_csv("../../annot/SRP065941/sample-annotation_stage1.tsv", separator="\t").filter(pl.col("sample_type") == "Tumor")



vcf_paths = sorted(glob.glob("../../vcf/SRP065941/vcf_filtered_pass_orientation/*/*.vcf"))
bam_paths = sorted(glob.glob(f"../../data/SRP065941/bam/*/*.bam"))


lookup_table = lookup_table.with_columns([
	pl.col("sample_name").map_elements(
		lambda id: [os.path.abspath(path) for path in vcf_paths if id in path][0],
		return_dtype=pl.String
	).alias("vcf_path"),
  
 	pl.col("sample_name").map_elements(
		lambda name : [os.path.abspath(path) for path in bam_paths if name in path][0],
		return_dtype=pl.String
	).alias("bam_path")
])

lookup_table_ffpe = lookup_table.filter(pl.col("preservation") == "FFPE")
# lookup_table_ffpe

# %%
for i in range(lookup_table_ffpe.shape[0]):
	sample_name = lookup_table_ffpe[i, "sample_name"]
	
	outdir = f"{outdir_root}/mut_info"
	os.makedirs(outdir, exist_ok=True)

	mut_info = make_mut_info(lookup_table_ffpe[i, "vcf_path"], sample_name, genome, ct_only = True)

	mut_info.write_csv(f"{outdir}/{sample_name}.microsec.mut-info.tsv", separator="\t")
	print(f"Finished creating Mut-Info for - {i + 1}. {sample_name}")

# %%
### Create Sample Info for MicroSEC

mut_info_paths = sorted(glob.glob(f"{outdir_root}/mut_info/*.microsec.mut-info.tsv"))
ref_path = os.path.abspath('../../data/ref/Homo_sapiens_assembly38.fasta')
mut_info_suffix = ".microsec.mut-info.tsv"

sample_info = pd.DataFrame()

for i in range(len(mut_info_paths)):
	
	sample_name = os.path.basename(mut_info_paths[i]).split(mut_info_suffix)[0]
	
	lookup = lookup_table_ffpe.filter(pl.col("sample_name").str.contains(sample_name))

	sample_info.loc[i, "sample_name"] = sample_name
	sample_info.loc[i, "mutation information tsv file"] = os.path.abspath(mut_info_paths[i])
	sample_info.loc[i, "BAM file"] = lookup[0, "bam_path"]
	
	sample_info.loc[i, "read length"] = str(100) # majoriy read length in bam files
	sample_info.loc[i, "adapter sequence read 1"] = np.nan
	sample_info.loc[i, "optional: adapter sequence read 2"] = np.nan
	sample_info.loc[i, "sample type: Human or Mouse"] = "hg38"
	sample_info.loc[i, "panel name"] = "TOP"
	sample_info.loc[i, "optional: reference genome fasta file"] = ref_path
	#sample_info.loc[i, "optional: simple repeat region bed file"] = np.nan

sample_info.to_csv(f"{outdir_root}/microsec.sample_info.tsv", sep="\t", index=False, header=False)


# %%
## PRJEB44073

print("Making inputs for ENA PRJEB44073")

outdir_root = "../PRJEB44073/vcf_filtered_pass_orientation/microsec/inputs"
genome = pysam.FastaFile("../../data/ref/Homo_sapiens_assembly38.fasta")
lookup_table = pl.read_csv("../../annot/PRJEB44073/sample-info_stage2.tsv", separator="\t")


vcf_paths = sorted(glob.glob("../../vcf/PRJEB44073/vcf_filtered_pass_orientation/*/*.vcf"))
bam_paths = sorted(glob.glob(f"../../data/PRJEB44073/bam/*/*.bam"))

# Skipping failed sample for now
lookup_table = lookup_table.filter(pl.col("sample_alias") != "Sample_B83_0029")

lookup_table = lookup_table.with_columns([
	pl.col("sample_alias").map_elements(
		lambda id: [os.path.abspath(path) for path in vcf_paths if id in path][0],
		return_dtype=pl.String
	).alias("vcf_path"),
  
 	pl.col("sample_alias").map_elements(
		lambda name : [os.path.abspath(path) for path in bam_paths if name in path][0],
		return_dtype=pl.String
	).alias("bam_path")
])

lookup_table_ffpe = lookup_table.filter(pl.col("preservation") == "FFPE")

# %%
for i in range(lookup_table_ffpe.shape[0]):
	sample_name = lookup_table_ffpe[i, "sample_alias"]
	
	outdir = f"{outdir_root}/mut_info"
	os.makedirs(outdir, exist_ok=True)

	mut_info = make_mut_info(lookup_table_ffpe[i, "vcf_path"], sample_name, genome, ct_only = True)

	mut_info.write_csv(f"{outdir}/{sample_name}.microsec.mut-info.tsv", separator="\t")
	print(f"Finished creating Mut-Info for - {i + 1}. {sample_name}")

# %%
### Create Sample Info for MicroSEC

mut_info_paths = sorted(glob.glob(f"{outdir_root}/mut_info/*.microsec.mut-info.tsv"))
ref_path = os.path.abspath('../../data/ref/Homo_sapiens_assembly38.fasta')
mut_info_suffix = ".microsec.mut-info.tsv"

sample_info = pd.DataFrame()

for i in range(len(mut_info_paths)):
	
	sample_name = os.path.basename(mut_info_paths[i]).split(mut_info_suffix)[0]
	
	lookup = lookup_table_ffpe.filter(pl.col("sample_alias").str.contains(sample_name))

	sample_info.loc[i, "sample_name"] = sample_name
	sample_info.loc[i, "mutation information tsv file"] = os.path.abspath(mut_info_paths[i])
	sample_info.loc[i, "BAM file"] = lookup[0, "bam_path"]
	
	sample_info.loc[i, "read length"] = str(98) # majoriy read length in bam files
	sample_info.loc[i, "adapter sequence read 1"] = np.nan
	sample_info.loc[i, "optional: adapter sequence read 2"] = np.nan
	sample_info.loc[i, "sample type: Human or Mouse"] = "hg38"
	sample_info.loc[i, "panel name"] = "TOP"
	sample_info.loc[i, "optional: reference genome fasta file"] = ref_path
	#sample_info.loc[i, "optional: simple repeat region bed file"] = np.nan

sample_info.to_csv(f"{outdir_root}/microsec.sample_info.tsv", sep="\t", index=False, header=False)



