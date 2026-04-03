#!/usr/bin/env Rscript
library(io)
library(precrec)
library(dplyr)
source("../common-ffpe-snvf/R/eval.R")

#######################################################

#' Evaluate FFPE Filter's Performance
#'
#' Evaluates the performance of an artifact filter on FFPE tumor samples by comparing
#' model's predictions against ground truth variants from matched frozen tumor samples.
#'
#' @param ffpe_tumoral Data frame of FFPE tumor samples with columns: case_id, sample_name
#' @param frozen_tumoral Data frame of frozen tumor samples with columns: case_id, vcf_fid
#' @param model_name Character string specifying the model/filter name
#' @param ff_vcf_dir Character string specifying the directory containing frozen tumor VCF files
#' @param ffpe_snvf_dir Character string specifying the directory containing FFPE filtering model's results
#' @param case_id_col Column name for case IDs (default: "case_id")
#' @param sample_id_col Column name for sample IDs (default: "sample_name")
#' @param snvf_ext File extension for score files (default: "tsv")
#'
#' @return Invisibly returns NULL. Writes per-sample and overall evaluation results to disk.
#'
#' @details
#' For each sample, the function:
#' 1. Retrieves FFPE and matched frozen tumor metadata
#' 2. Reads model predictions and ground truth variants
#' 3. Preprocesses the data using model-specific logic
#' 4. Evaluates model performance if both true and false labels exist
#' 5. Writes results to disk
#'
#' Samples are skipped if no overlap or complete overlap exists between FFPE and frozen variants.
evaluate_sample_set <- function(
	ffpe_tumoral,
	frozen_tumoral,
	model_name,
	ffpe_snvf_dir,
	case_id_col = "case_id",
	sample_id_col = "sample_name",
	snvf_ext = "snv"
) {

	dataset <- basename(dirname(ffpe_snvf_dir))
	variant_set <- basename(ffpe_snvf_dir)
	message(sprintf("Evaluating %s on Dataset: %s | Variant Set: %s:", toupper(model_name), dataset, variant_set))
	outdir_root <- file.path(dataset, variant_set)

	for (i in seq_len(nrow(ffpe_tumoral))){

		## Get FFPE sample metadata
		case_id <- ffpe_tumoral[i, case_id_col]
		sample_name <- ffpe_tumoral[i, sample_id_col]

		message(sprintf("	%d. Processing Sample: %s", i, sample_name))

		snvf_path <- file.path(ffpe_snvf_dir, model_name, sample_name, sprintf("%s.%s.%s", sample_name, model_name, snvf_ext))
		if (!file.exists(snvf_path)){
			message(sprintf("		Warning: %s SNVF does not exist at %s . Skipping", model_name, snvf_path))
			next
		}

		## Read prepared ground truth
		ground_truth_path <- file.path("..", "ground-truth", dataset, gsub("-micr1234", "", variant_set), sample_name, sprintf("%s.ground-truth.tsv", sample_name))
		ground_truth <- read.delim(ground_truth_path)
		ground_truth$truth <- as.logical(ground_truth$truth)

		# Read model's score for current sample and Apply model specific processing
		d <- preprocess_filter(read.delim(snvf_path), model_name)
		d <- merge(ground_truth, d, by=c("chrom", "pos", "ref", "alt"))

		## Ensure that the ground truth annotation are not all true or all false
		## These if so they cannot be evaluated
		## Write the scores with truth labels to disk for overall evaluation and skip evaluation
		if((nrow(d[d$truth, ]) == 0)){
			message(sprintf("		No true labels exist for %s. Skipping individual evaluation...", snvf_path))
			write_sample_eval(d, NULL, outdir_root, sample_name, model_name)
			next
		}
		if((nrow(d[!d$truth, ]) == 0)){
			message(sprintf("		No false labels exist for %s. Skipping individual evaluation...", snvf_path))
			write_sample_eval(d, NULL, outdir_root, sample_name, model_name)
			next
		}

		# Evaluate the filter's performance
		res <- evaluate_filter(d, model_name, sample_name)

		# Write results
		write_sample_eval(d, res, outdir_root, sample_name, model_name)
		
	}

	# Overall Evaluation
	## The scores annotated with ground truth is combined into a single dataframe
	message("	Performing overall evaluation across all samples")

	all_score_truth <- bind_rows(
		lapply(ffpe_tumoral[[sample_id_col]], function(sample_name) {
			path <- file.path(outdir_root, "model-scores_truths", sample_name, sprintf("%s_%s-scores_truths.tsv", sample_name, model_name))
			if (!file.exists(path)){
				message(sprintf("		Warning: %s was not was not found. SKIPPING", path))
			} else {
				d <- read.delim(path)
				d$alt <- as.character(d$alt)
				d$sample_name <- sample_name
				return(d)
			}
		})
	)

	# Evaluate across all samples
	overall_res <- evaluate_filter(all_score_truth, model_name, "all-samples")
	write_overall_eval(all_score_truth, overall_res, outdir_root, "all-samples", model_name)
	message(cat("\tComplete.\n"))

}


########################################################
# Evaluate the ffpe filter
model_name <- "vafsnvf"

# #########################################  ENA PRJEB8754  #######################################

# # Read Annotation Table
# lookup_table <- read.delim("../annot/PRJEB8754/sample-info_matched-ff-ffpe_on-pat-id-sample-type.tsv")

# # Stratify annotation table based on FFPE and FF Somatic Variants
# ffpe_tumoral <- lookup_table[(lookup_table$preservation == "FFPE"), ]
# frozen_tumoral <- lookup_table[(lookup_table$preservation == "Frozen"), ]


# ## Evaluate dataset
# evaluate_sample_set(
# 	ffpe_tumoral = ffpe_tumoral,
# 	frozen_tumoral = frozen_tumoral,
# 	model_name = model_name,
# 	ffpe_snvf_dir = "../ffpe-snvf/PRJEB8754/dup-unmarked_raw_filtermutectcalls_obmm_unfiltered"
# )

# evaluate_sample_set(
# 	ffpe_tumoral = ffpe_tumoral,
# 	frozen_tumoral = frozen_tumoral,
# 	model_name = model_name,
# 	ffpe_snvf_dir = "../ffpe-snvf/PRJEB8754/dup-unmarked_filtered_pass-orient-pos-sb"
# )

# evaluate_sample_set(
# 	ffpe_tumoral = ffpe_tumoral,
# 	frozen_tumoral = frozen_tumoral,
# 	model_name = model_name,
# 	ffpe_snvf_dir = "../ffpe-snvf/PRJEB8754/dup-unmarked_filtered_pass-orient-pos-sb-vaf-dp"
# )

# evaluate_sample_set(
# 	ffpe_tumoral = ffpe_tumoral,
# 	frozen_tumoral = frozen_tumoral,
# 	model_name = model_name,
# 	ffpe_snvf_dir = "../ffpe-snvf/PRJEB8754/dup-unmarked_filtered_pass-orient-pos-sb-vaf-dp-macni"
# )

# # #################################  ENA PRJEB44073  ########################################

# # Read Annotation Table
# lookup_table <- read.delim("../annot/PRJEB44073/sample-info_stage2.tsv")
# lookup_table$sample_name <- lookup_table$sample_alias

# # Stratify annotation table based on FFPE and FF Somatic Variants
# ffpe_tumoral <- lookup_table[(lookup_table$preservation == "FFPE"), ]
# frozen_tumoral <- lookup_table[(lookup_table$preservation == "Frozen"), ]

# ## Evaluate  dataset
# evaluate_sample_set(
# 	ffpe_tumoral = ffpe_tumoral,
# 	frozen_tumoral = frozen_tumoral,
# 	model_name = model_name,
# 	ffpe_snvf_dir = "../ffpe-snvf/PRJEB44073/filtered_pass-orientation"
# )


# evaluate_sample_set(
# 	ffpe_tumoral = ffpe_tumoral,
# 	frozen_tumoral = frozen_tumoral,
# 	model_name = model_name,
# 	ffpe_snvf_dir = "../ffpe-snvf/PRJEB44073/filtered_pass-orientation-exome"
# )

# evaluate_sample_set(
# 	ffpe_tumoral = ffpe_tumoral,
# 	frozen_tumoral = frozen_tumoral,
# 	model_name = model_name,
# 	ffpe_snvf_dir = "../ffpe-snvf/PRJEB44073/filtered_pass-orientation-exome-blacklist"
# )

# evaluate_sample_set(
# 	ffpe_tumoral = ffpe_tumoral,
# 	frozen_tumoral = frozen_tumoral,
# 	model_name = model_name,
# 	ffpe_snvf_dir = "../ffpe-snvf/PRJEB44073/filtered_pass-orientation-exome-blacklist-macni"
# )

# evaluate_sample_set(
# 	ffpe_tumoral = ffpe_tumoral,
# 	frozen_tumoral = frozen_tumoral,
# 	model_name = model_name,
# 	ffpe_snvf_dir = "../ffpe-snvf/PRJEB44073/filtered_pass-orientation-exome-blacklist-macni-micr1234"
# )


# # ##################################  ENA SRP044740  ##############################################

# Read Annotation Table
lookup_table <- read.delim("../annot/SRP044740/sample-info_stage2.tsv")

# Stratify annotation table based on FFPE and FF Somatic Variants
ffpe_tumoral <- lookup_table[(lookup_table$preservation == "FFPE"), ]
frozen_tumoral <- lookup_table[(lookup_table$preservation == "Frozen"), ]

# ## Evaluate  dataset
# evaluate_sample_set(
# 	ffpe_tumoral = ffpe_tumoral,
# 	frozen_tumoral = frozen_tumoral,
# 	model_name = model_name,
# 	ffpe_snvf_dir = "../ffpe-snvf/SRP044740/filtered_pass-orientation"
# )


evaluate_sample_set(
	ffpe_tumoral = ffpe_tumoral,
	frozen_tumoral = frozen_tumoral,
	model_name = model_name,
	ffpe_snvf_dir = "../ffpe-snvf/SRP044740/filtered_pass-orientation-exome"
)

evaluate_sample_set(
	ffpe_tumoral = ffpe_tumoral,
	frozen_tumoral = frozen_tumoral,
	model_name = model_name,
	ffpe_snvf_dir = "../ffpe-snvf/SRP044740/filtered_pass-orientation-exome-blacklist"
)

evaluate_sample_set(
	ffpe_tumoral = ffpe_tumoral,
	frozen_tumoral = frozen_tumoral,
	model_name = model_name,
	ffpe_snvf_dir = "../ffpe-snvf/SRP044740/filtered_pass-orientation-exome-blacklist-macni"
)

evaluate_sample_set(
	ffpe_tumoral = ffpe_tumoral,
	frozen_tumoral = frozen_tumoral,
	model_name = model_name,
	ffpe_snvf_dir = "../ffpe-snvf/SRP044740/filtered_pass-orientation-exome-blacklist-macni-micr1234"
)


##################################  ENA SRP065941  ##############################################

# Read Annotation Table
lookup_table <- read.delim("../annot/SRP065941/sample-annotation_stage1.tsv")

# Stratify annotation table based on FFPE and FF Somatic Variants
ffpe_tumoral <- lookup_table[(lookup_table$preservation == "FFPE") & (lookup_table$sample_type == "Tumor"), ]
frozen_tumoral <- lookup_table[(lookup_table$preservation == "frozen")  & (lookup_table$sample_type == "Tumor"), ]


# ## Evaluate dataset
# evaluate_sample_set(
# 	ffpe_tumoral = ffpe_tumoral,
# 	frozen_tumoral = frozen_tumoral,
# 	model_name = model_name,
# 	ffpe_snvf_dir = "../ffpe-snvf/SRP065941/filtered_pass-orientation"
# )

evaluate_sample_set(
	ffpe_tumoral = ffpe_tumoral,
	frozen_tumoral = frozen_tumoral,
	model_name = model_name,
	ffpe_snvf_dir = "../ffpe-snvf/SRP065941/filtered_pass-orientation-exome"
)

evaluate_sample_set(
	ffpe_tumoral = ffpe_tumoral,
	frozen_tumoral = frozen_tumoral,
	model_name = model_name,
	ffpe_snvf_dir = "../ffpe-snvf/SRP065941/filtered_pass-orientation-exome-blacklist"
)

evaluate_sample_set(
	ffpe_tumoral = ffpe_tumoral,
	frozen_tumoral = frozen_tumoral,
	model_name = model_name,
	ffpe_snvf_dir = "../ffpe-snvf/SRP065941/filtered_pass-orientation-exome-blacklist-micr1234"
)
