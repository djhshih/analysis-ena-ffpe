#!/usr/bin/env Rscript
library(io)
library(precrec)
source("../common-ffpe-snvf/R/eval.R")

#######################################################

#' Load or Build Truth Set for a Case and Variant Caller
#'
#' This function loads a previously cached truth set for a given case and variant caller,
#' or constructs a new truth set if one does not already exist. The truth set is created
#' by taking the union of variants from matched fresh-frozen (FF) VCF files and cached by
#' assigning unique identifiers to each variant.
#'
#' @param case_id Character. The identifier for the case/tumor being analyzed.
#' @param variant_caller Character. The name of the variant caller used to generate
#'   the VCF files (e.g., "mutect2", "varscan2").
#' @param matched_ff_vcf_paths Character vector. File paths to the matched fresh-frozen
#'   VCF files to be used for creating the truth set.
#' @param outdir Character. The output directory where truth sets will be cached.
#'   Defaults to \code{main.outdir}.
#'
#' @return A data frame or tibble containing the truth set with unique variant identifiers.
#'   If cached, returns the previously saved truth set. Otherwise, returns a newly
#'   constructed truth set created from the union of variants in matched_ff_vcf_paths.
#'
#' @details
#' Truth sets are cached in the \code{truth_sets} subdirectory of \code{outdir}
#' with the naming convention: \code{<case_id>_<variant_caller>.tsv}
#'
#' @examples
#' \dontrun{
#' truth <- get_truth_set("CASE001", "mutect2", c("file1.vcf", "file2.vcf"))
get_truth_set <- function(case_id, matched_ff_vcf_paths, outdir = main.outdir) {
	truth_set_dir <- file.path(outdir, "truth_sets")
	dir.create(truth_set_dir, showWarnings = FALSE, recursive = TRUE)
	truth_set_path <- file.path(
		truth_set_dir, 
		paste(sapply(strsplit(basename(matched_ff_vcf_paths), "\\."), `[`, 1), collapse = "_"), 
		sprintf("%s.tsv", case_id)
	)
	
	if (file.exists(truth_set_path)){
		message(cat("\tGround truth set exists, reading from file"))
		truth <- qread(truth_set_path)
	} else {
		message(cat("\tCached ground truth set does not exists, generating ground truth set from:"))
		message(cat("\t", matched_ff_vcf_paths))
		truth <- snv_union(matched_ff_vcf_paths)
		truth <- add_id(truth)
		qwrite(truth, truth_set_path)
	}
	return(truth)
}

#' Evaluate a Set of Tumor Samples
#'
#' This function evaluates the performance of a variant calling model on a set of FFPE (Formalin-Fixed Paraffin-Embedded)
#' tumor samples by comparing model predictions against a ground truth set derived from matched frozen tumor samples.
#' It processes each sample individually, then performs an overall evaluation across all samples.
#'
#' @param ffpe_tumoral A data frame containing FFPE tumor sample metadata. Must include columns:
#'   - case_id: character, unique identifier for the case
#'   - tumor_bam_fid: character, file ID for the tumor BAM file
#'   - workflow_type: character, name of the variant caller used
#'
#' @param frozen_tumoral A data frame containing frozen tumor sample metadata. Must include columns:
#'   - case_id: character, unique identifier matching cases in ffpe_tumoral
#'   - vcf_fid: character, file ID for the VCF file
#'
#' @param model_name A character string specifying the name of the model being evaluated.
#'   Used for file naming and directory organization.
#'
#' @param outdir_root A character string specifying the root output directory where evaluation
#'   results will be written.
#'
#' @param ff_vcf_dir A character string specifying the directory containing frozen tumor VCF files.
#'
#' @param ffpe_snvf_dir A character string specifying the directory containing FFPE model score files.
#'   Expected structure: ffpe_snvf_dir/model_name/case_id/
#'
#' @return Invisibly returns NULL. The function writes evaluation results to disk:
#'   - Per-sample evaluation results
#'   - Overall evaluation across all samples
#'
#' @details
#' For each sample, the function:
#' 1. Retrieves FFPE and matched frozen tumor metadata
#' 2. Reads model predictions and ground truth variants
#' 3. Preprocesses the data using model-specific logic
#' 4. Evaluates model performance if both true and false labels exist
#' 5. Writes results to disk
#'
#' Samples are skipped if no overlap exists between FFPE and frozen variants or if
#' either true or false truth labels are absent.
evaluate_sample_set <- function(
	ffpe_tumoral,
	frozen_tumoral,
	model_name,
	ff_vcf_dir,
	ffpe_snvf_dir,
	case_id_col = "case_id",
	sample_id_col = "sample_name"
) {

	for (i in seq_len(nrow(ffpe_tumoral))){

		## Get FFPE sample metadata
		case_id <- ffpe_tumoral[i, case_id_col]
		sample_name <- ffpe_tumoral[i, sample_id_col]
		variant_caller <- "Mutect2"
		dataset <- basename(dirname(ffpe_snvf_dir))
		variant_set <- basename(ffpe_snvf_dir)
		outdir_root <- file.path(dataset, variant_set)

		snvf_path <- file.path(ffpe_snvf_dir, model_name, sample_name, sprintf("%s.%s.tsv", sample_name, model_name))
		if (!file.exists(snvf_path)){
			message(sprintf("	Warning: %s SNVF does not exist at %s . Skipping", model_name, snvf_path))
			next
		}

		## Read prepared ground truth
		truth <- read.delim(file.path("..", "ground-truth", dataset, gsub("-micr1234", "", variant_set), sample_name, sprintf("%s.ground-truth.tsv", sample_name)))
		truth$truth <- as.logical(truth$truth)

		# read model's score for current sample and Apply model specific processing
		d <- preprocess_ffpolish(read.delim(snvf_path))
		d <- merge(d, truth, by=c("chrom", "pos", "ref", "alt"))

		## Check if true labels exist in the variant_score_truth table (d)
		## If not this means there's no overlap between FFPE and FF variants
		## Cases like these are skipped as evaluation is not supported by precrec
		if((nrow(d[d$truth, ]) == 0)){
			message(sprintf("	no true labels exist for %s", snvf_path))
			write_sample_eval(d, NULL, outdir_root, sample_name, model_name)
			next
		}
		if((nrow(d[!d$truth, ]) == 0)){
			message(sprintf("	no false labels exist for %s", snvf_path))
			write_sample_eval(d, NULL, outdir_root, sample_name, model_name)
			next
		}

		message(sprintf("	%s", snvf_path))

		# Evaluate the filter's performance
		res <- evaluate_filter(d, model_name, sample_name)

		# Write results
		write_sample_eval(d, res, outdir_root, sample_name, model_name)
		
	}

	# Overall Evaluation
	## The scores annotated with ground truth is combined into a single dataframe
	message("	performing Evaluation across all samples")

	all_score_truth <- do.call(
		rbind,
		lapply(ffpe_tumoral[[sample_id_col]], function(sample_name) {
			path <- file.path(outdir_root, "model-scores_truths", sample_name, sprintf("%s_%s-scores_truths.tsv", sample_name, model_name))
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
	overall_res <- evaluate_filter(all_score_truth, model_name, "all-samples")
	write_overall_eval(all_score_truth, overall_res, outdir_root, "all-samples", model_name)

}


########################################################
# Evaluate the ffpe filter
model_name <- "ffpolish"
message(sprintf("Evaluating %s:", model_name))


# #########################################  ENA PRJEB8754  #######################################

# Read Annotation Table
lookup_table <- read.delim("../annot/PRJEB8754/sample-info_matched-ff-ffpe_on-pat-id-sample-type.tsv")

# Stratify annotation table based on FFPE and FF Somatic Variants
ffpe_tumoral <- lookup_table[(lookup_table$preservation == "FFPE"), ]
frozen_tumoral <- lookup_table[(lookup_table$preservation == "Frozen"), ]


# ## Evaluate dataset
# evaluate_sample_set(
# 	ffpe_tumoral = ffpe_tumoral,
# 	frozen_tumoral = frozen_tumoral,
# 	model_name = model_name,
# 	ff_vcf_dir = "../vcf/PRJEB8754/raw_filtermutectcalls_obmm_unfiltered_dup-unmarked",
# 	ffpe_snvf_dir = "../ffpe-snvf/PRJEB8754/raw_filtermutectcalls_obmm_unfiltered_dup-unmarked"
# )


# ## Evaluate dataset
# evaluate_sample_set(
# 	ffpe_tumoral = ffpe_tumoral,
# 	frozen_tumoral = frozen_tumoral,
# 	model_name = model_name,
# 	ff_vcf_dir = "../vcf/PRJEB8754/filtered_pass-orient-pos-sb-ad-blacklist-macni_dup-unmarked",
# 	ffpe_snvf_dir = "../ffpe-snvf/PRJEB8754/filtered_pass-orient-pos-sb-ad-blacklist-macni_dup-unmarked"
# )

#################################  ENA PRJEB44073  ########################################

# Read Annotation Table
lookup_table <- read.delim("../annot/PRJEB44073/sample-info_stage2.tsv")
lookup_table$sample_name <- lookup_table$sample_alias

# Stratify annotation table based on FFPE and FF Somatic Variants
ffpe_tumoral <- lookup_table[(lookup_table$preservation == "FFPE"), ]
frozen_tumoral <- lookup_table[(lookup_table$preservation == "Frozen"), ]

## Evaluate  dataset
evaluate_sample_set(
	ffpe_tumoral = ffpe_tumoral,
	frozen_tumoral = frozen_tumoral,
	model_name = model_name,
	ff_vcf_dir = "../vcf/PRJEB44073/filtered_pass-orientation-dp20",
	ffpe_snvf_dir = "../ffpe-snvf/PRJEB44073/filtered_pass-orientation-dp20"
)

## Evaluate  dataset
evaluate_sample_set(
	ffpe_tumoral = ffpe_tumoral,
	frozen_tumoral = frozen_tumoral,
	model_name = model_name,
	ff_vcf_dir = "../vcf/PRJEB44073/filtered_pass-orientation-dp20-blacklist",
	ffpe_snvf_dir = "../ffpe-snvf/PRJEB44073/filtered_pass-orientation-dp20-blacklist"
)

## Evaluate  dataset
evaluate_sample_set(
	ffpe_tumoral = ffpe_tumoral,
	frozen_tumoral = frozen_tumoral,
	model_name = model_name,
	ff_vcf_dir = "../vcf/PRJEB44073/filtered_pass-orientation-dp20-blacklist-macni",
	ffpe_snvf_dir = "../ffpe-snvf/PRJEB44073/filtered_pass-orientation-dp20-blacklist-macni"
)

## Evaluate  dataset
evaluate_sample_set(
	ffpe_tumoral = ffpe_tumoral,
	frozen_tumoral = frozen_tumoral,
	model_name = model_name,
	ff_vcf_dir = "../vcf/PRJEB44073/filtered_pass-orientation-dp20-blacklist-macni",
	ffpe_snvf_dir = "../ffpe-snvf/PRJEB44073/filtered_pass-orientation-dp20-blacklist-macni-micr1234"
)


##################################  ENA SRP044740  ##############################################

# Read Annotation Table
lookup_table <- read.delim("../annot/SRP044740/sample-info_stage2.tsv")

# Stratify annotation table based on FFPE and FF Somatic Variants
ffpe_tumoral <- lookup_table[(lookup_table$preservation == "FFPE"), ]
frozen_tumoral <- lookup_table[(lookup_table$preservation == "Frozen"), ]

## Evaluate  dataset
evaluate_sample_set(
	ffpe_tumoral = ffpe_tumoral,
	frozen_tumoral = frozen_tumoral,
	model_name = model_name,
	ff_vcf_dir = "../vcf/SRP044740/filtered_pass-orientation-dp20",
	ffpe_snvf_dir = "../ffpe-snvf/SRP044740/filtered_pass-orientation-dp20"
)

## Evaluate  dataset
evaluate_sample_set(
	ffpe_tumoral = ffpe_tumoral,
	frozen_tumoral = frozen_tumoral,
	model_name = model_name,
	ff_vcf_dir = "../vcf/SRP044740/filtered_pass-orientation-dp20-blacklist",
	ffpe_snvf_dir = "../ffpe-snvf/SRP044740/filtered_pass-orientation-dp20-blacklist"
)

## Evaluate  dataset
evaluate_sample_set(
	ffpe_tumoral = ffpe_tumoral,
	frozen_tumoral = frozen_tumoral,
	model_name = model_name,
	ff_vcf_dir = "../vcf/SRP044740/filtered_pass-orientation-dp20-blacklist-macni",
	ffpe_snvf_dir = "../ffpe-snvf/SRP044740/filtered_pass-orientation-dp20-blacklist-macni"
)

## Evaluate  dataset
evaluate_sample_set(
	ffpe_tumoral = ffpe_tumoral,
	frozen_tumoral = frozen_tumoral,
	model_name = model_name,
	ff_vcf_dir = "../vcf/SRP044740/filtered_pass-orientation-dp20-blacklist-macni",
	ffpe_snvf_dir = "../ffpe-snvf/SRP044740/filtered_pass-orientation-dp20-blacklist-macni-micr1234"
)


##################################  ENA SRP065941  ##############################################

# Read Annotation Table
lookup_table <- read.delim("../annot/SRP065941/sample-annotation_stage1.tsv")

# Stratify annotation table based on FFPE and FF Somatic Variants
ffpe_tumoral <- lookup_table[(lookup_table$preservation == "FFPE") & (lookup_table$sample_type == "Tumor"), ]
frozen_tumoral <- lookup_table[(lookup_table$preservation == "frozen")  & (lookup_table$sample_type == "Tumor"), ]


## Evaluate dataset
evaluate_sample_set(
	ffpe_tumoral = ffpe_tumoral,
	frozen_tumoral = frozen_tumoral,
	model_name = model_name,
	ff_vcf_dir = "../vcf/SRP065941/filtered_pass-orientation-dp20",
	ffpe_snvf_dir = "../ffpe-snvf/SRP065941/filtered_pass-orientation-dp20"
)


evaluate_sample_set(
	ffpe_tumoral = ffpe_tumoral,
	frozen_tumoral = frozen_tumoral,
	model_name = model_name,
	ff_vcf_dir = "../vcf/SRP065941/filtered_pass-orientation-dp20-blacklist",
	ffpe_snvf_dir = "../ffpe-snvf/SRP065941/filtered_pass-orientation-dp20-blacklist"
)

evaluate_sample_set(
	ffpe_tumoral = ffpe_tumoral,
	frozen_tumoral = frozen_tumoral,
	model_name = model_name,
	ff_vcf_dir = "../vcf/SRP065941/filtered_pass-orientation-dp20-blacklist",
	ffpe_snvf_dir = "../ffpe-snvf/SRP065941/filtered_pass-orientation-dp20-blacklist-micr1234"
)
