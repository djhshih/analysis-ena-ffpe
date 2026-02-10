#!/usr/bin/env Rscript
library(io)
source("../../common-ffpe-snvf/R/macni_somatic.R")

## Filtering Function
process_dataset <- function(dataset, variant_set, new_variant_set){

	message("Processing dataset: ", dataset)
	vcf_paths = Sys.glob(sprintf("../../vcf/%s/%s/*/*.vcf", dataset, variant_set))

	for (i in seq_along(vcf_paths)){
		path <- vcf_paths[i]
		sample_name <- basename(dirname(path))
		message(sprintf("	%s. sample: %s", i, sample_name))

		res <- run_macni(path, sample_name, thresh = 0.5)
		somatic <- res[res$is_somatic, c("chrom", "pos", "ref", "alt")]

		outdir <- file.path("..", dataset, new_variant_set, sample_name)
		dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

		write.table(somatic, file.path(outdir, sprintf("%s.tsv", sample_name)), sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
	}
}

## Filter Datasets
process_dataset(
	dataset = "PRJEB8754",
	variant_set = "vcf_filtered_pass-orient-pos-sb-ad-blacklist_dup-unmarked",
	new_variant_set = "filtered_pass-orient-pos-sb-ad-blacklist-macni_dup-unmarked"
)

process_dataset(
	dataset = "PRJEB44073",
	variant_set = "vcf_filtered_pass-orientation-dp10-blacklist",
	new_variant_set = "filtered_pass-orientation-dp10-blacklist-macni"
)

process_dataset(
	dataset = "SRP044740",
	variant_set = "vcf_filtered_pass-orientation-dp10-blacklist",
	new_variant_set = "filtered_pass-orientation-dp10-blacklist-macni"
)


