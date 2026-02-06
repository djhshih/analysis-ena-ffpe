#!/usr/bin/env Rscript
library(io)
library(ggplot2)
library(arrow)
source("../../common-ffpe-snvf/R/macni_somatic.R")

#' Make per-sample VAF plots
#'
#' @description
#' Create histogram plots of variant allele frequency (or another numeric column)
#' facetted by sample. Plots are filled by the `is_somatic` flag and use free y-scales
#' per sample. The plot may be returned as a ggplot object or written to a file via qdraw.
#'
#' @param all_res (data.frame) Data frame containing at minimum: `sample_name`, `is_somatic`,
#'   and the column specified by `x_col` (default "vaf").
#' @param x_col (string) Name of the column to plot on the x-axis. Default: "vaf".
#' @param title (string|null) Plot title.
#' @param subtitle (string|null) Plot subtitle.
#' @param x_lab (string) Label for x-axis. Default: "Variant Allele Frequency".
#' @param y_lab (string) Label for y-axis. Default: "count".
#' @param file (string|null) If provided, path to save the plot using qdraw(). If NULL, the ggplot object is returned.
#' @param base_size (numeric) Base font size passed to theme_classic(). Default: 20.
#' @param bins (integer) Number of bins for geom_histogram(). Default: 50.
#'
#' @return (ggplot object) If `file` is NULL the ggplot object is returned. If `file` is provided,
#'   the plot is saved to disk and nothing is returned.
#'
make_vaf_plot <- function(
		all_res, 
		x_col = "vaf", 
		title = NULL, 
		subtitle = NULL, 
		x_lab = "Variant Allele Frequency", 
		y_lab = "count",
		file = NULL,
		base_size = 20, 
		bins = 50,
		width = NULL,
		height = NULL
	){

	if(is.null(width) & is.null(height)){
		options(repr.plot.width = floor(length(unique(all_res$sample_name))/1.5), repr.plot.height = 20)
	} else {
		options(repr.plot.width = width, repr.plot.height = height)
	}

	macni_vaf_plot <- ggplot(na.omit(all_res), aes(x = .data[[x_col]], fill = is_somatic)) +
		geom_histogram(bins = bins, color="black") +
		facet_wrap(~sample_name, ncol = 4, scales = "free_y") +
		theme_classic(base_size = base_size) +
		labs(
			title = title,
			subtitle = subtitle,
			x = x_lab,
			y = y_lab
		)

	if(!is.null(file)){
		qdraw(macni_vaf_plot, file, width = 20, height = floor(length(unique(all_res$sample_name))/1.5))
	} else {
		return(macni_vaf_plot)
	}
}

## ENA PRJEB44073 | Chong21
vcf_paths = Sys.glob("../../vcf/PRJEB44073/vcf_filtered_pass-orientation-dp10-blacklist/*/*.vcf")
all_res <- list()

for (i in seq_along(vcf_paths)){
	path <- vcf_paths[i]
	sample_name <- basename(dirname(path))
	message(sprintf("%s. Processing sample: %s", i, sample_name))

	vcf <- qread(path)
	res <- run_macni(vcf,  sample_name, thresh = 0.5)
	res$sample_name <- sample_name
	
	all_res[[i]] <- res
}

all_res <- do.call(rbind, all_res)
write_parquet(all_res, file.path("ena_chong21_macni_ppcut-0.5_table.parquet"))

vcf$Sample_B83_0001[[1]]$F1R2

thresholds <- seq(0, 1, 0.1)

outdir = "vaf_plots"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

for (thresh in thresholds){
	message("Making VAF plots using posterior probability cutoff theshold: ", thresh)

	all_res$is_somatic <- ifelse(all_res$macni_pp > thresh, TRUE, FALSE)

	make_vaf_plot(
		all_res, 
		title = "ENA Chong21", 
		subtitle = sprintf("MACNI Posterior Probability threshold > %s", thresh),
		file = file.path(outdir, sprintf("ena_chong21_macni-%s_per-sample-vaf-plots.pdf", thresh))
	)
}

## ENA SRP044740 
vcf_paths = Sys.glob("../../vcf/SRP044740/vcf_filtered_pass-orientation-dp10-blacklist/*/*.vcf")
all_res <- list()

for (i in seq_along(vcf_paths)){
	path <- vcf_paths[i]
	sample_name <- basename(dirname(path))
	message(sprintf("Processing sample: %s", sample_name))

	vcf <- qread(path)
	res <- run_macni(vcf,  sample_name, thresh = 0.5)
	res$sample_name <- sample_name
	
	all_res[[i]] <- res
}

all_res <- do.call(rbind, all_res)
write_parquet(all_res, "ena_SRP044740_macni_ppcut-0.5_table.parquet")


thresholds <- seq(0, 1, 0.1)

outdir = "vaf_plots"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

for (thresh in thresholds){
	message("Making VAF plots using posterior probability cutoff theshold: ", thresh)

	all_res$is_somatic <- ifelse(all_res$macni_pp > thresh, TRUE, FALSE)

	make_vaf_plot(
		all_res, 
		title = "ENA Chong21", 
		subtitle = sprintf("MACNI Posterior Probability threshold > %s", thresh),
		file = file.path(outdir, sprintf("ena_SRP044740_macni-%s_per-sample-vaf-plots.pdf", thresh))
	)
}

