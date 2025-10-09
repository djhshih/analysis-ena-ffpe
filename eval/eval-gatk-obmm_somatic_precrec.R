library(io)
library(precrec)

# TODO move common functions into module
# TODO refactoring in progress

source("../common-ffpe-snvf/R/eval.R")

#######################################################

## Function to lookup matched Fresh Frozen sample, read snvs and annotate id
# @params lookup_table (data.frame) | sample annotation table
# @params vcf_dir (string) | Root directory where the vcfs are located
# @return data.frame | Table of annotated ground truth snv
get_ground_truth <- function(matched_ff_sample_name, vcf_dir) {
    matched_ff <- read_vcf(file.path(vcf_dir, matched_ff_sample_name, sprintf("%s.vcf", matched_ff_sample_name)), columns = c("chrom", "pos", "ref", "alt"))
    matched_ff <- add_id(matched_ff)
    matched_ff
}

model_name <- "gatk-obmm"

########################################################
message("Evaluating gatk_obmm:")

# Setup Directories
## specific for each dataset

dataset_id <- "PRJEB8754"
message(sprintf("Dataset: %s", dataset_id))
# Directory for FFPE SNVF inputs
ffpe_snvf.dir <- sprintf("../ffpe-snvf/%s/vcf_pass-orient-pos-sb_ad_filtered", dataset_id)
# Directory for the Somatic VCFs
vcf.dir <- sprintf("../vcf/%s/vcf_pass-orient-pos-sb_ad_filtered", dataset_id)
# Output directory
main.outdir <- sprintf("%s/vcf_pass-orient-pos-sb_ad_filtered", dataset_id)
score_truth_outdir <- sprintf("%s/vcf_pass-orient-pos-sb_ad_filtered/model-scores_truths", dataset_id)
eval_outdir <- sprintf("%s/vcf_pass-orient-pos-sb_ad_filtered/roc-prc-auc/precrec", dataset_id)


# Read Annotation Table
lookup_table <- read.delim(sprintf("../annot/%s/sample-info_matched-ff-ffpe_on-pat-id-sample-type.tsv", dataset_id))

# Stratify annotation table based on FFPE and FF Somatic Variants
ffpe_tumoral <- lookup_table[(lookup_table$preservation == "FFPE"), ]
frozen_tumoral <- lookup_table[(lookup_table$preservation == "Frozen"), ]


# Evaluate gatk_obmm
## Per Sample
## Each sample is evaluated first due to the necessity of independently annotating the scores with ground truth

message("Processing samples:")
for (index in seq_len(nrow(ffpe_tumoral))){

	## Get FFPE sample metadata
	sample_name <- sprintf("%s_%s",ffpe_tumoral[index, "sample_alias"], ffpe_tumoral[index, "run_accession"])
	pat_id <- ffpe_tumoral[index, "inferred_id"]

	## Getting matched FF metadata by matching patient ID
	matched_ff_metadata <- frozen_tumoral[(frozen_tumoral$inferred_id == pat_id), ]
    matched_ff_sample_name <- sprintf("%s_%s",matched_ff_metadata[1, "sample_alias"], matched_ff_metadata[1, "run_accession"])

	d <- read_gatk_snv(sample_name, model_name, vcf.dir)
	truth <- get_ground_truth(matched_ff_sample_name, vcf.dir)
	d <- preprocess_gatk_obmm(d, truth)

	## Check if true labels exist in the variant_score_truth table (d)
	## If not this means there's no overlap between FFPE and FF variants
	## Cases like these are skipped as evaluation is not supported by precrec
	if(nrow(d[d$truth, ]) == 0){
		message(sprintf("	no true labels exist for %s", sample_name))
		next
	}

	message(sprintf("	%s", sample_name))

	# Evaluate the filter's performance
	gatk_obmm_res <- evaluate_filter(d, model_name, sample_name)

	# write results
	write_sample_eval(d, gatk_obmm_res, main.outdir, sample_name, model_name)

}



# Overall Evaluation
## The scores annotated with ground truth is combined into a single dataframe
message("	performing Evaluation across all samples")
gatk_obmm_all_score_truth <- do.call(
	rbind,
	lapply(seq_len(nrow(ffpe_tumoral)), function(i) {
		sample_name <- sprintf("%s_%s",ffpe_tumoral[i, "sample_alias"], ffpe_tumoral[i, "run_accession"])
		path <- file.path(score_truth_outdir, sample_name, sprintf("%s_%s-scores_truths.tsv", sample_name, model_name))
		if (!file.exists(path)){
			message(sprintf("	Warning: %s was not was not found. SKIPPING", path))
		} else {
			d <- read.delim(path)
			d$sample_name <- sample_name
			d
		}
	})
)


# Evaluate across all samples
gatk_obmm_overall_res <- evaluate_filter(gatk_obmm_all_score_truth, model_name, "all-samples")
write_overall_eval(gatk_obmm_all_score_truth, gatk_obmm_overall_res, score_truth_outdir, eval_outdir, "all-samples", model_name)


#############################################

# Setup Directories and lookup table for ENA SRP044740

dataset_id <- "SRP044740"
message(sprintf("Dataset: %s", dataset_id))
# Directory for FFPE SNVF inputs
ffpe_snvf.dir <- sprintf("../ffpe-snvf/%s/vcf_filtered_pass_orientation", dataset_id)
# Directory for the Somatic VCFs
vcf.dir <- sprintf("../vcf/%s/vcf_filtered_pass_orientation", dataset_id)
# Output directory
main.outdir <- sprintf("%s/vcf_filtered_pass_orientation", dataset_id)
score_truth_outdir <- sprintf("%s/vcf_filtered_pass_orientation/model-scores_truths", dataset_id)
eval_outdir <- sprintf("%s/vcf_filtered_pass_orientation/roc-prc-auc/precrec", dataset_id)


# Read Annotation Table
lookup_table <- read.delim(sprintf("../annot/%s/sample-info_stage1.tsv", dataset_id))

# Stratify annotation table based on FFPE and FF Somatic Variants
ffpe_tumoral <- lookup_table[(lookup_table$sample_type == "FFPE"), ]
frozen_tumoral <- lookup_table[(lookup_table$sample_type == "FROZ"), ]

# Evaluate gatk_obmm
## Per Sample
## Each sample is evaluated first due to the necessity of independently annotating the scores with ground truth
message("Processing samples:")
for (index in seq_len(nrow(ffpe_tumoral))){

	## Get FFPE sample metadata
	sample_name <- sprintf("%s_%s",ffpe_tumoral[index, "sample_alias"], ffpe_tumoral[index, "run_accession"])
	pat_id <- ffpe_tumoral[index, "sample_number"]
	message(sprintf("	%s", sample_name))

	## Getting matched FF metadata by matching patient ID
	matched_ff_metadata <- frozen_tumoral[(frozen_tumoral$sample_number == pat_id), ]
    matched_ff_sample_name <- sprintf("%s_%s",matched_ff_metadata[1, "sample_alias"], matched_ff_metadata[1, "run_accession"])

	d <- read_gatk_snv(sample_name, model_name, vcf.dir)
	truth <- get_ground_truth(matched_ff_sample_name, vcf.dir)
	d <- preprocess_gatk_obmm(d, truth)

	# Evaluate the filter's performance
	gatk_obmm_res <- evaluate_filter(d, model_name, sample_name)

	# write results
	write_sample_eval(d, gatk_obmm_res, main.outdir, sample_name, model_name)

}

# Overall Evaluation
## The scores annotated with ground truth is combined into a single dataframe
message("	performing Evaluation across all samples")
gatk_obmm_all_score_truth <- do.call(
	rbind,
	lapply(seq_len(nrow(ffpe_tumoral)), function(i) {
		sample_name <- sprintf("%s_%s",ffpe_tumoral[i, "sample_alias"], ffpe_tumoral[i, "run_accession"])
		path <- file.path(score_truth_outdir, sample_name, sprintf("%s_%s-scores_truths.tsv", sample_name, model_name))
		if (!file.exists(path)){
			message(sprintf("	Warning: %s was not was not found. SKIPPING", path))
		} else {
			d <- read.delim(path)
			d$sample_name <- sample_name
			d
		}
	})
)



# Evaluate across all samples
gatk_obmm_overall_res <- evaluate_filter(gatk_obmm_all_score_truth, model_name, "all-samples")
write_overall_eval(gatk_obmm_all_score_truth, gatk_obmm_overall_res, score_truth_outdir, eval_outdir, "all-samples", model_name)


############################################

# Setup Directories and lookup table for ENA SRP065941

dataset_id <- "SRP065941"
message(sprintf("Dataset: %s", dataset_id))
# Directory for FFPE SNVF inputs
ffpe_snvf.dir <- sprintf("../ffpe-snvf/%s/vcf_filtered_pass_orientation", dataset_id)
# Directory for the Somatic VCFs
vcf.dir <- sprintf("../vcf/%s/vcf_filtered_pass_orientation", dataset_id)
# Output directory
main.outdir <- sprintf("%s/vcf_filtered_pass_orientation", dataset_id)
score_truth_outdir <- sprintf("%s/vcf_filtered_pass_orientation/model-scores_truths", dataset_id)
eval_outdir <- sprintf("%s/vcf_filtered_pass_orientation/roc-prc-auc/precrec", dataset_id)


# Read Annotation Table
lookup_table <- read.delim(sprintf("../annot/%s/sample-annotation_stage1.tsv", dataset_id))

# Stratify annotation table based on FFPE and FF Somatic Variants
ffpe_tumoral <- lookup_table[(lookup_table$preservation == "FFPE") & (lookup_table$sample_type == "Tumor"), ]
frozen_tumoral <- lookup_table[(lookup_table$preservation == "frozen") & (lookup_table$sample_type == "Tumor"), ]

# Evaluate gatk_obmm
## Per Sample
## Each sample is evaluated first due to the necessity of independently annotating the scores with ground truth
message("Processing samples:")
for (index in seq_len(nrow(ffpe_tumoral))){

	## Get FFPE sample metadata
	sample_name <- ffpe_tumoral[index, "sample_name"]
	case_id <- ffpe_tumoral[index, "case_id"]
	message(sprintf("	%s", sample_name))

	## Getting matched FF metadata by matching patient ID
	matched_ff_metadata <- frozen_tumoral[(frozen_tumoral$case_id == case_id), ]
    matched_ff_sample_name <- matched_ff_metadata[1, "sample_name"]

	d <- read_gatk_snv(sample_name, model_name, vcf.dir)
	truth <- get_ground_truth(matched_ff_sample_name, vcf.dir)
	d <- preprocess_gatk_obmm(d, truth)

	# Evaluate the filter's performance
	gatk_obmm_res <- evaluate_filter(d, model_name, sample_name)

	# write results
	write_sample_eval(d, gatk_obmm_res, main.outdir, sample_name, model_name)

}

# Overall Evaluation
## The scores annotated with ground truth is combined into a single dataframe
message("	performing Evaluation across all samples")
gatk_obmm_all_score_truth <- do.call(
	rbind,
	lapply(seq_len(nrow(ffpe_tumoral)), function(i) {
		sample_name <- ffpe_tumoral[i, "sample_name"]
		path <- file.path(score_truth_outdir, sample_name, sprintf("%s_%s-scores_truths.tsv", sample_name, model_name))
		if (!file.exists(path)){
			message(sprintf("	Warning: %s was not was not found. SKIPPING", path))
		} else {
			d <- read.delim(path)
			d$sample_name <- sample_name
			d
		}
	})
)



# Evaluate across all samples
gatk_obmm_overall_res <- evaluate_filter(gatk_obmm_all_score_truth, model_name, "all-samples")
write_overall_eval(gatk_obmm_all_score_truth, gatk_obmm_overall_res, score_truth_outdir, eval_outdir, "all-samples", model_name)



#############################################

# Setup Directories and lookup table for ENA SRP065941

dataset_id <- "PRJEB44073"
message(sprintf("Dataset: %s", dataset_id))
# Directory for FFPE SNVF inputs
ffpe_snvf.dir <- sprintf("../ffpe-snvf/%s/vcf_filtered_pass_orientation", dataset_id)
# Directory for the Somatic VCFs
vcf.dir <- sprintf("../vcf/%s/vcf_filtered_pass_orientation", dataset_id)
# Output directory
main.outdir <- sprintf("%s/vcf_filtered_pass_orientation", dataset_id)
score_truth_outdir <- sprintf("%s/vcf_filtered_pass_orientation/model-scores_truths", dataset_id)
eval_outdir <- sprintf("%s/vcf_filtered_pass_orientation/roc-prc-auc/precrec", dataset_id)


# Read Annotation Table
lookup_table <- read.delim(sprintf("../annot/%s/sample-info_stage2.tsv", dataset_id))
lookup_table <- lookup_table[lookup_table$sample_alias != "Sample_B83_0029", ]

# Stratify annotation table based on FFPE and FF Somatic Variants
ffpe_tumoral <- lookup_table[(lookup_table$preservation == "FFPE"), ]
frozen_tumoral <- lookup_table[(lookup_table$preservation == "Frozen"), ]

# Evaluate gatk_obmm
## Per Sample
## Each sample is evaluated first due to the necessity of independently annotating the scores with ground truth
message("Processing samples:")
for (index in seq_len(nrow(ffpe_tumoral))){

	## Get FFPE sample metadata
	sample_name <- ffpe_tumoral[index, "sample_alias"]
	case_id <- ffpe_tumoral[index, "case_id"]
	message(sprintf("	%s", sample_name))

	## Getting matched FF metadata by matching patient ID
	matched_ff_metadata <- frozen_tumoral[(frozen_tumoral$case_id == case_id), ]
    matched_ff_sample_name <- matched_ff_metadata[1, "sample_alias"]

	d <- read_gatk_snv(sample_name, model_name, vcf.dir)
	truth <- get_ground_truth(matched_ff_sample_name, vcf.dir)
	d <- preprocess_gatk_obmm(d, truth)

	# Evaluate the filter's performance
	gatk_obmm_res <- evaluate_filter(d, model_name, sample_name)

	# write results
	write_sample_eval(d, gatk_obmm_res, main.outdir, sample_name, model_name)

}

# Overall Evaluation
## The scores annotated with ground truth is combined into a single dataframe
message("	performing Evaluation across all samples")
gatk_obmm_all_score_truth <- do.call(
	rbind,
	lapply(seq_len(nrow(ffpe_tumoral)), function(i) {
		sample_name <- ffpe_tumoral[i, "sample_alias"]
		path <- file.path(score_truth_outdir, sample_name, sprintf("%s_%s-scores_truths.tsv", sample_name, model_name))
		if (!file.exists(path)){
			message(sprintf("	Warning: %s was not was not found. SKIPPING", path))
		} else {
			d <- read.delim(path)
			d$sample_name <- sample_name
			d
		}
	})
)



# Evaluate across all samples
gatk_obmm_overall_res <- evaluate_filter(gatk_obmm_all_score_truth, model_name, "all-samples")
write_overall_eval(gatk_obmm_all_score_truth, gatk_obmm_overall_res, score_truth_outdir, eval_outdir, "all-samples", model_name)

