import os
import glob
import polars as pl

## Functions
def return_path_if_exists(path: str, abs: bool=False) -> str:
	if os.path.exists(path):
		return os.path.abspath(path) if abs else path
	else:
		raise FileNotFoundError(f"Path {path} does not exist.")

def create_filtering_scripts(models: list, vcf_path: str, bam_path: str, filtered_outdir: str, sample_name: str, ref_path: str) -> None:
	for model in models:
	
		if model == "ideafix":
			
			outdir = "script_ideafix"
			os.makedirs(outdir, exist_ok=True)
					
			## Template script for each model
			template_dir = return_path_if_exists(f"../../common-ffpe-snvf/R/{model}.R", abs=True)
			content = f"#!/usr/bin/env Rscript\nRscript {template_dir} --vcf '{vcf_path}' --ref '{ref_path}' --outdir '{filtered_outdir}'\n"
			filename = f"{outdir}/{model}_{sample_name}.sh"
			
			with open(filename, "w") as f:
				f.write(content)

		elif model == "ffpolish":
			outdir = "script_ffpolish"
			os.makedirs(outdir, exist_ok=True)

			content = [
				"#!/usr/bin/env bash \n",
				f"ffpolish filter -o {filtered_outdir}/ffpolish/{sample_name} -p {sample_name} {ref_path} {vcf_path} {bam_path} \n"
			]
			filename = f"{outdir}/{model}_{sample_name}.sh"
			
			with open(filename, "w") as f:
				f.writelines(content)
				
		else:
			
			outdir = "script"
			os.makedirs(outdir, exist_ok=True)
			
			template_dir = return_path_if_exists(f"../../common-ffpe-snvf/templates/ffpe-snvf/{model}.sh.template", abs=True)
		
			content = f"#!/bin/bash\nbash {template_dir} '{bam_path}' '{vcf_path}' '{filtered_outdir}/{model}/{sample_name}' '{ref_path}'\n"
			filename = f"{outdir}/{model}_{sample_name}.sh"
		
			with open(filename, "w") as f:
				f.write(content)
			
		print(f"Created {model} filtering script for {vcf_path}")

#--------------------

## user parameters
models = ["mobsnvf", "vafsnvf", "sobdetector", "ideafix", "ffpolish"]
ref_path = return_path_if_exists("../../data/ref/Homo_sapiens_assembly38.fasta", abs=True)

# Process each dataset

## PRJEB8754
vcf_paths = sorted([os.path.abspath(path) for path in glob.glob("../../vcf/PRJEB8754/raw_filtermutectcalls_obmm_unfiltered_dup-unmarked/*/*.vcf")])
filtered_outdir = os.path.abspath("../PRJEB8754/raw_filtermutectcalls_obmm_unfiltered_dup-unmarked")
bam_dir = "../../data/PRJEB8754/bam_dup-unmarked"

# skip non FFPE vcfs
for vcf_path in vcf_paths:
	if "Frozen" in vcf_path:
		continue
	
	sample_name = os.path.basename(vcf_path).split(".")[0]
	
	bam_path = return_path_if_exists(f"{bam_dir}/{sample_name}/{sample_name}.bam", abs=True)
	create_filtering_scripts(models, vcf_path, bam_path, filtered_outdir, sample_name, ref_path)


## PRJEB44073
vcf_paths = sorted([os.path.abspath(path) for path in glob.glob("../../vcf/PRJEB44073/filtered_pass-orientation/*/*.vcf")])
filtered_outdir = os.path.abspath("../PRJEB44073/filtered_pass-orientation")
bam_dir = "../../data/PRJEB44073/bam"

annot = pl.read_csv("../../annot/PRJEB44073/sample-info_stage3.tsv", separator = "\t")
ffpe_samples = annot.filter(pl.col("preservation") == "FFPE").get_column("sample_name").to_list()

# skip non FFPE vcfs
vcf_paths = [path for path in vcf_paths if any(sample in path for sample in ffpe_samples)]

for vcf_path in vcf_paths:
	
	vcf_filename = os.path.basename(vcf_path)
	sample_name = vcf_filename.removesuffix(".vcf")
	
	bam_path = return_path_if_exists(f"{bam_dir}/{sample_name}/{sample_name}.bam", abs=True)
	create_filtering_scripts(models, vcf_path, bam_path, filtered_outdir, sample_name, ref_path)
	
print("All filtering scripts created.")


## SRP044740
vcf_paths = sorted([os.path.abspath(path) for path in glob.glob("../../vcf/SRP044740/filtered_pass-orientation/*/*.vcf")])
filtered_outdir = os.path.abspath("../SRP044740/filtered_pass-orientation")
bam_dir = "../../data/SRP044740/bam"

# skip non FFPE vcfs
for vcf_path in vcf_paths:
	if "FROZ" in vcf_path:
		continue
	
	vcf_filename = os.path.basename(vcf_path)
	sample_name = vcf_filename.removesuffix(".vcf")
	
	bam_path = return_path_if_exists(f"{bam_dir}/{sample_name}.bam", abs=True)
	create_filtering_scripts(models, vcf_path, bam_path, filtered_outdir, sample_name, ref_path)


## SRP065941
vcf_paths = sorted([os.path.abspath(path) for path in glob.glob("../../vcf/SRP065941/filtered_pass-orientation/*/*.vcf")])
filtered_outdir = os.path.abspath("../SRP065941/filtered_pass-orientation")
bam_dir = "../../data/SRP065941/bam"

# skip non FFPE vcfs
for vcf_path in vcf_paths:
	if "frozen" in vcf_path.lower():
		continue
	
	vcf_filename = os.path.basename(vcf_path)
	sample_name = vcf_filename.removesuffix(".vcf")
	
	bam_path = return_path_if_exists(f"{bam_dir}/{sample_name}/{sample_name}.bam", abs=True)
	create_filtering_scripts(models, vcf_path, bam_path, filtered_outdir, sample_name, ref_path)

