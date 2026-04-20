#!/usr/bin/env python
import polars as pl
import os
import glob
import sys

# Local Dependencies
sys.path.append("../../common-ffpe-snvf/python")
from microsec_preparation import *
from common import return_path_if_exists

### Main Execution
#### Resource path and Directory Setup
outdir_root = "."
repo_root = ".."
vcf_root = "../vcf"
annot_root = f"{repo_root}/annot"
batch_script_root = "script_microsec"

simple_repeat_path = return_path_if_exists(f"{repo_root}/data/ref/ucsc_simple-repeat_hg38.bed", abs=True)
ref_path = return_path_if_exists(f"{repo_root}/data/ref/Homo_sapiens_assembly38.fasta", abs=True)
microsec_script_path = return_path_if_exists(f"{repo_root}/common-ffpe-snvf/R/microsec.R", abs=True)

## Wrapper function to prepare inputs and execution script for a dataset
def prepare_dataset_inputs(
		dataset: str,
		variant_set: str,
		bam_vcf_paths: pl.DataFrame,
		vcf_path_col: str = "vcf_path", 
		bam_path_col: str = "bam_path",
		sample_name_col: str = "sample_name",
		simple_repeat_path: str = simple_repeat_path,
		microsec_script_path: str = microsec_script_path,
		ref_path: str = ref_path,
		batch_script_root: str = batch_script_root,
		ct_only: bool = False,
		snv_only: bool = False,
		lazy: bool = False
	):

	microsec_root = f"{outdir_root}/{dataset}/{variant_set}"

	print(f"Creating MicroSEC input for Dataset: {dataset}, Variant Set: {variant_set}")

	for i in range(bam_vcf_paths.height):

		## sample specific variable setup
		sample_name = bam_vcf_paths[i, sample_name_col]
		vcf_path = bam_vcf_paths[i, vcf_path_col]
		bam_path = bam_vcf_paths[i, bam_path_col]

		print(f"\t{i+1}. Processing Sample: {sample_name}, from VCF: {vcf_path}")

		###--------------------------------------------


		# Read in the simple repeats bed file to annotate variants later
		simple_repeats_bed = pl.read_csv(
			simple_repeat_path, 
			separator="\t", 
			has_header=False, 
			new_columns=["chrom", "start", "end", "type", "score"],
		).select(["chrom", "start", "end"])
		
		###--------------------------------------------
		
		## Prepare mut_info for the sample
		output_dir = f"{microsec_root}/{sample_name}/inputs"
		os.makedirs(output_dir, exist_ok=True)
		mut_info_path = os.path.abspath(f"{output_dir}/{sample_name}.microsec.mut-info.tsv")

		if not (lazy and os.path.exists(mut_info_path)):
			print("\t\t- Preparing MicroSEC Mutation Info...")
			mut_info = make_mut_info(
				vcf_path=vcf_path, 
				ref_fasta_path=ref_path, 
				sample_name=sample_name,
				simple_repeats_bed=simple_repeats_bed,
				ct_only=ct_only, 
				snv_only=snv_only
			)
		
			## Write mut_info
			mut_info.write_csv(mut_info_path, separator="\t")
			print(f"\t\t\tPrepared MicroSEC Mutation-Info with {mut_info.height} variants")
			print(f"\t\t\tMicroSEC Mutation-Info saved to {mut_info_path}")

		###----------------------------------------------

		sample_info_path = os.path.abspath(f"{microsec_root}/{sample_name}/inputs/{sample_name}.microsec.sample_info.tsv")
		
		## Prepare sample_info for the sample
		if not (lazy and os.path.exists(sample_info_path)):
			print("\t\t- Preparing sample info")
			sample_info = make_sample_info(sample_name, mut_info_path, bam_path, ref_path, simple_repeat_path)

			## Write sample_info
			pl.DataFrame(sample_info).write_csv(sample_info_path, separator="\t")

		###-----------------------------------------------

		# Write batch execution scripts
		print("\t\t- Preparing execution script")
		batch_script_dir = f"{batch_script_root}/{dataset}"
		os.makedirs(batch_script_dir, exist_ok=True)
		script_path = f"{batch_script_dir}/{sample_name}.microsec.sh"
		microsec_outdir = os.path.abspath(f"{microsec_root}/{sample_name}")
		
		with open(script_path, "w") as file:
			file.writelines("#!/usr/bin/env bash\n")
			file.writelines(f"Rscript {microsec_script_path} --sample_info '{sample_info_path}' --output_dir '{microsec_outdir}'\n")

	print("Done.\n")



# Prepare microsec inputs and batch execution scripts

# #### Prepare PRJEB44073
# ## Create annotation table linking ffpe vcfs to bams
# dataset = "PRJEB44073"
# variant_set = "filtered_pass-orientation-exome"
# ffpe_vcf_paths = [os.path.abspath(path) for path in sorted(glob.glob(f"{vcf_root}/{dataset}/{variant_set}/*/*.vcf.gz"))]

# bam_vcf_paths_df = (
# 	pl.DataFrame({
# 		"sample_name" : [path.split("/")[-2] for path in ffpe_vcf_paths],
# 		"vcf_path" : ffpe_vcf_paths,
# 		"bam_path" : [return_path_if_exists(f"{repo_root}/data/{dataset}/bam/{path.split("/")[-2]}/{path.split("/")[-2]}.bam", abs=True) for path in ffpe_vcf_paths],
# 	})
# 	.sort("sample_name")
# )

# prepare_dataset_inputs(
# 	dataset = dataset, 
# 	variant_set = variant_set,
# 	bam_vcf_paths = bam_vcf_paths_df
# )


#### Prepare SRP044740
## Create annotation table linking ffpe samples to 
dataset = "SRP044740"
variant_set = "filtered_pass-orientation-exome-blacklist-macni"
ffpe_vcf_paths = [os.path.abspath(path) for path in sorted(glob.glob(f"{vcf_root}/{dataset}/{variant_set}/*/*.vcf.gz"))]

sample_paths = pl.DataFrame({
		"sample_name" : [path.split("/")[-2] for path in ffpe_vcf_paths],
		"vcf_path" : ffpe_vcf_paths,
		"bam_path" : [return_path_if_exists(f"{repo_root}/data/{dataset}/bam/{path.split("/")[-2]}.bam", abs=True) for path in ffpe_vcf_paths],
	})

## Prepare microsec inputs and batch execution scripts
prepare_dataset_inputs(
	dataset = f"{dataset}", 
	variant_set = variant_set,
	bam_vcf_paths = sample_paths,
	ct_only=True
)


# #### Prepare SRP065941
# ## Create annotation table linking ffpe samples to 
# dataset = "SRP065941"
# variant_set = "filtered_pass-orientation"
# ffpe_vcf_paths = [os.path.abspath(path) for path in sorted(glob.glob(f"{vcf_root}/{dataset}/{variant_set}/*/*.vcf.gz"))]

# sample_paths = pl.DataFrame({
# 		"sample_name" : [path.split("/")[-2] for path in ffpe_vcf_paths],
# 		"vcf_path" : ffpe_vcf_paths,
# 		"bam_path" : [return_path_if_exists(f"{repo_root}/data/{dataset}/bam/{path.split("/")[-2]}/{path.split("/")[-2]}.bam", abs=True) for path in ffpe_vcf_paths],
# 	})

# ## Prepare microsec inputs and batch execution scripts
# prepare_dataset_inputs(
# 	dataset = f"{dataset}", 
# 	variant_set = variant_set, 
# 	bam_vcf_paths = sample_paths,
# 	ct_only=True
# )
