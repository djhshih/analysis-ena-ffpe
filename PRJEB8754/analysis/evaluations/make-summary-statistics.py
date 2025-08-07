#!/usr/bin/python3
import polars as pl
import os
import glob

root_outdir = "../../evaluations/vcf_pass-orient-pos-sb_ad_filtered"
vcf_dir = "../../data/vcf_pass-orient-pos-sb_ad_filtered"
ffpe_snvf_dir="../../ffpe-snvf/vcf_pass-orient-pos-sb_ad_filtered"

def import_formatted_vcf(vcf_path, sample_name, snv_only = True, standard_chr_only=True) -> pl.DataFrame:
	
	"""
	Import a VCF or SNV file and return a Polars DataFrame with standardized columns.

	Parameters
	----------
	vcf_path : str
		Path to the VCF or SNV file.
	snv_only : bool, optional
		If True, filter to include only single nucleotide variants (default: True).

	Returns
	-------
	pl.DataFrame
		DataFrame with columns: chrom, pos, ref, alt, sample_id.
		Multi-allelic variants are exploded into separate rows.
	"""
	
	allowed_chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
	
	# print(sample_name)
	
	if vcf_path.endswith(".vcf") or vcf_path.endswith(".vcf.gz"):
		vcf = (
			pl.read_csv(vcf_path, separator="\t", comment_prefix="##", columns=["#CHROM", "POS", "REF", "ALT"])
			# .filter(pl.col("#CHROM").is_in(allowed_chroms))
			.rename({
						"#CHROM" : "chrom",
						"POS" : "pos", 
						"REF": "ref", 
						"ALT" : "alt"
			})
		)
	elif vcf_path.endswith(".snv"):
		vcf = pl.read_csv(vcf_path, separator="\t", infer_schema_length=10000)
		vcf = vcf.rename({col: col.lower() for col in vcf.columns})
		
	
	vcf = (
		vcf.with_columns([
			pl.lit(sample_name).alias("sample_id"),
			pl.col("ref").str.to_uppercase(),
			pl.col("alt").str.to_uppercase()
		])
		# .filter(pl.col("chrom").is_in(allowed_chroms))
		.with_columns(
			pl.col("alt").str.split(",")
		)
		.explode("alt")
	)
	
	if standard_chr_only:
		vcf = vcf.filter(pl.col("chrom").is_in(allowed_chroms))
 
	if snv_only:
		vcf = vcf.filter((pl.col("alt").str.len_chars() == 1) & (pl.col("ref").str.len_chars() == 1))
	
	return vcf

c_to_t_mask = (
	((pl.col("ref") == "C") & (pl.col("alt") == "T")) | 
	((pl.col("ref") == "G") & (pl.col("alt") == "A"))
)

snv_mask = (
	(pl.col("alt").str.len_chars() == 1) & 
	(pl.col("ref").str.len_chars() == 1)
)

## List paths to VCFs
all_vcf_paths = glob.glob(f"{vcf_dir}/*/*.vcf")

## Read in the sample annotation table
lookup_table = pl.read_csv("../../annot/sample-info_matched-ff-ffpe_on-pat-id.tsv", separator="\t")

## Make a dataframe with only the FFPE samples
ffpe_vcf_ids = lookup_table.filter(pl.col("preservation") == "FFPE").get_column("sample_alias")

## Filter the list to only keep the FFPE VCF paths
ffpe_vcf_paths = [vcf_path for vcf_path in all_vcf_paths if any(fid in vcf_path for fid in ffpe_vcf_ids)]

# Get mutation counts for the FFPE vcfs and groud truth. [SNVs, Non-SNVs, C>T Mutations, Non-C>T mutations]
metrics = []

for path in ffpe_vcf_paths:
	sample_id = path.split("/")[-2]
	patient_id = sample_id.split("_")[0]
	sample_alias = "_".join(sample_id.split("_")[:3])
	run_acc = sample_id.split("_")[-1]

	truth_annot = (
     	lookup_table
    	.filter(
         	(pl.col("inferred_id") == patient_id) &
			(pl.col("sample_type") == "Meta") &
			(pl.col("preservation") == "Frozen")
        ))

	
	truth_sample_id = f"{truth_annot[0, "sample_alias"]}_{truth_annot[0, "run_accession"]}"
	variant_caller = "MuTect2"
	
	ffpe_vcf = import_formatted_vcf(path, sample_alias, snv_only=False)
	total_mutations = ffpe_vcf.shape[0]
	total_snvs = ffpe_vcf.filter(snv_mask).shape[0]
	ct_snvs = ffpe_vcf.filter(snv_mask).filter(c_to_t_mask).shape[0]
	non_ct_snvs = ffpe_vcf.filter(snv_mask).filter(~c_to_t_mask).shape[0]


	ground_truth = import_formatted_vcf(f"{vcf_dir}/{truth_sample_id}/{truth_sample_id}.vcf", sample_name = truth_sample_id)
	
	gt_total_mutations = ground_truth.shape[0]
	gt_total_snvs = ground_truth.filter(snv_mask).shape[0]
	gt_ct_snvs = ground_truth.filter(snv_mask).filter(c_to_t_mask).shape[0]
	gt_non_ct_snvs = ground_truth.filter(snv_mask).filter(~c_to_t_mask).shape[0]

	ffpe_snvs_in_gt = (
		ffpe_vcf
	 	.filter(snv_mask)
	   	.join(ground_truth, on=["chrom", "pos", "ref" , "alt"], how = "semi")
	).shape[0]
 
	ffpe_snvs_not_in_gt = (
		ffpe_vcf
	 	.filter(snv_mask)
	   	.join(ground_truth, on=["chrom", "pos", "ref" , "alt"], how = "anti")
	).shape[0]

	ffpe_ct_snvs_in_gt = (
		ffpe_vcf
	 	.filter(snv_mask)
	  	.filter(c_to_t_mask)
	   	.join(ground_truth, on=["chrom", "pos", "ref" , "alt"], how = "semi")
	).shape[0]
 
	ffpe_ct_snvs_not_in_gt = (
		ffpe_vcf
	 	.filter(snv_mask)
	  	.filter(c_to_t_mask)
	   	.join(ground_truth, on=["chrom", "pos", "ref" , "alt"], how = "anti")
	).shape[0]

	ffpe_nct_snvs_in_gt = (
		ffpe_vcf
	 	.filter(snv_mask)
	  	.filter(~c_to_t_mask)
	   	.join(ground_truth, on=["chrom", "pos", "ref" , "alt"], how = "semi")
	).shape[0]
 
	ffpe_nct_snvs_not_in_gt = (
		ffpe_vcf
	 	.filter(snv_mask)
	  	.filter(~c_to_t_mask)
	   	.join(ground_truth, on=["chrom", "pos", "ref" , "alt"], how = "anti")
	).shape[0]

	# Ground truth SNVs in FFPE and not in FFPE
	gt_snvs_in_ffpe = (
		ground_truth
		.filter(snv_mask)
		.join(ffpe_vcf, on=["chrom", "pos", "ref", "alt"], how="semi")
	).shape[0]

	gt_snvs_not_in_ffpe = (
		ground_truth
		.filter(snv_mask)
		.join(ffpe_vcf, on=["chrom", "pos", "ref", "alt"], how="anti")
	).shape[0]

	# Ground truth C>T SNVs in FFPE and not in FFPE
	gt_ct_snvs_in_ffpe = (
		ground_truth
		.filter(snv_mask)
		.filter(c_to_t_mask)
		.join(ffpe_vcf, on=["chrom", "pos", "ref", "alt"], how="semi")
	).shape[0]

	gt_ct_snvs_not_in_ffpe = (
		ground_truth
		.filter(snv_mask)
		.filter(c_to_t_mask)
		.join(ffpe_vcf, on=["chrom", "pos", "ref", "alt"], how="anti")
	).shape[0]

	# Ground truth non-C>T SNVs in FFPE and not in FFPE
	gt_nct_snvs_in_ffpe = (
		ground_truth
		.filter(snv_mask)
		.filter(~c_to_t_mask)
		.join(ffpe_vcf, on=["chrom", "pos", "ref", "alt"], how="semi")
	).shape[0]

	gt_nct_snvs_not_in_ffpe = (
		ground_truth
		.filter(snv_mask)
		.filter(~c_to_t_mask)
		.join(ffpe_vcf, on=["chrom", "pos", "ref", "alt"], how="anti")
	).shape[0]

	metrics.append({
		"ffpe_sample_id": sample_id,
  		"variant_caller": variant_caller,
		"ffpe_total_mutations": total_mutations,
		"ffpe_total_non_snvs": total_mutations - total_snvs,
		"ffpe_total_snvs": total_snvs,
		"ffpe_ct_snvs": ct_snvs,
		"ffpe_non_ct_snvs": non_ct_snvs,
		"gt_total_mutations": gt_total_mutations,
		"gt_total_non_snvs": gt_total_mutations - gt_total_snvs,
		"gt_total_snvs": gt_total_snvs,
		"gt_ct_snvs": gt_ct_snvs,
		"gt_non_ct_snvs": gt_non_ct_snvs,
		"ffpe_snvs_in_gt" : ffpe_snvs_in_gt,
		"ffpe_snvs_not_in_gt" : ffpe_snvs_not_in_gt,
		"ffpe_ct_snvs_in_gt" : ffpe_ct_snvs_in_gt,
		"ffpe_ct_snvs_not_in_gt" : ffpe_ct_snvs_not_in_gt,
		"ffpe_non_ct_snvs_in_gt" : ffpe_nct_snvs_in_gt,
		"ffpe_non_ct_snvs_not_in_gt" : ffpe_nct_snvs_not_in_gt,
		"gt_snvs_in_ffpe": gt_snvs_in_ffpe,
		"gt_snvs_not_in_ffpe": gt_snvs_not_in_ffpe,
		"gt_ct_snvs_in_ffpe": gt_ct_snvs_in_ffpe,
		"gt_ct_snvs_not_in_ffpe": gt_ct_snvs_not_in_ffpe,
		"gt_non_ct_snvs_in_ffpe": gt_nct_snvs_in_ffpe,
		"gt_non_ct_snvs_not_in_ffpe": gt_nct_snvs_not_in_ffpe,
	})

metrics_df = pl.DataFrame(metrics).sort(["ffpe_sample_id"])

metrics_df.write_csv(f"{root_outdir}/ffpe_frozen_variants-summary_statistics.tsv", separator="\t")
print(f"Saved summary statistics for FFPE and Frozen (Ground Truth). \nLocation: {root_outdir}/ffpe_frozen_variants-summary_statistics.tsv\n")

### Make a metrics dataframe for outputs of each ffpe filtering model
def variant_count(variants) -> tuple:
	
	all_mutations = variants.shape[0]
	snvs_only = variants.filter(snv_mask).shape[0]
	non_snvs = variants.filter(~snv_mask).shape[0]
	c_to_t = variants.filter(c_to_t_mask).shape[0]
	non_c_to_t = variants.filter(~c_to_t_mask).shape[0]

	return all_mutations, snvs_only, non_snvs, c_to_t, non_c_to_t

counts = []

for path in ffpe_vcf_paths:
    
	sample_id = path.split("/")[-2]


	ffpe_vcf = import_formatted_vcf(path, sample_id, snv_only=False)
	total_mutations = ffpe_vcf.shape[0]
	total_snvs = ffpe_vcf.filter(snv_mask).shape[0]
	ct_snvs = ffpe_vcf.filter(snv_mask).filter(c_to_t_mask).shape[0]
	non_ct_snvs = ffpe_vcf.filter(snv_mask).filter(~c_to_t_mask).shape[0]

	mobsnvf_vcf = import_formatted_vcf(f"{ffpe_snvf_dir}/mobsnvf/{sample_id}/{sample_id}.mobsnvf.snv", sample_id, snv_only=False, standard_chr_only=True)
	vafsnvf_vcf = import_formatted_vcf(f"{ffpe_snvf_dir}/vafsnvf/{sample_id}/{sample_id}.vafsnvf.snv", sample_id, snv_only=False, standard_chr_only=True)
	sobdetector_vcf = import_formatted_vcf(f"{ffpe_snvf_dir}/sobdetector/{sample_id}/{sample_id}.sobdetector.snv", sample_id, snv_only=False, standard_chr_only=True)

	mobsnvf_mutations, mobsnvf_snvs, mobsnvf_non_snvs, mobsnvf_c_to_t, mobsnvf_non_c_to_t =  variant_count(mobsnvf_vcf)
	vafsnvf_mutations, vafsnvf_snvs, vafsnvf_non_snvs, vafsnvf_c_to_t, vafsnvf_non_c_to_t =  variant_count(vafsnvf_vcf)
	sobdetector_mutations, sobdetector_snvs, sobdetector_non_snvs, sobdetector_c_to_t, sobdetector_non_c_to_t =  variant_count(mobsnvf_vcf)


	counts.append({
		"ffpe_sample_id": sample_id,
  		"variant_caller": variant_caller,
		"ffpe_vcf_snvs": total_snvs,
		"mobsnvf_snvs" : mobsnvf_snvs,
		"vafsnvf_snvs" : vafsnvf_snvs,
		"sobdetector_snvs" : sobdetector_snvs,
		"mobsnvf_c_to_t_snvs" : mobsnvf_c_to_t,
		"ffpe_vcf_c_to_t_snvs": ct_snvs,
		"vafsnvf_c_to_t_snvs" : vafsnvf_c_to_t,
		"sobdetector_c_to_t_snvs" : sobdetector_c_to_t,
		"ffpe_vcf_non_ct_snvs": non_ct_snvs,
		"mobsnvf_non_c_to_t_snvf" : mobsnvf_non_c_to_t,
		"vafsnvf_non_c_to_t_snvf" : vafsnvf_non_c_to_t,
		"sobdetector_non_c_to_t_snvf" : sobdetector_non_c_to_t,
	})

model_variant_count = pl.DataFrame(counts).sort(["ffpe_sample_id"])

model_variant_count.write_csv(f"{root_outdir}/models_variants-summary_statistics.tsv", separator="\t")
print(f"Saved summary statistics for variants present across the outputs for each model \nLocation: {root_outdir}/models_variants-summary_statistics.tsv\n")
