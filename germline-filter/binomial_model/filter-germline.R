#!/usr/bin/env Rscript
library(io)
library(dplyr)
library(tidyr)
library(ggplot2)
library(parallel)
library(arrow)
options(mc.cores=4)


## Functions
binom_tests <- function(xs, ns, p, alternative) {
	unlist(mcmapply(
		function(x, n, p, alternative) {
			if(any(is.na(c(x, n)))){
				NA
			} else {
				binom.test(x, n, p, alternative)$p.value
			}
		},
		xs, ns, p, alternative,
		SIMPLIFY=FALSE
	))
}


## Germline Filtering
dataset_variants <- c(
	"SRP044740/SRP044740_all-variants_filtered_pass-orientation-dp10-blacklist.parquet",
	"PRJEB44073/PRJEB44073_all-variants_filtered_pass-orientation-dp10-blacklist.parquet"

)

for (variant_set in dataset_variants) {

	message(sprintf("Processing variants from %s", variant_set))
	dataset <- basename(dirname(variant_set))
	filename <- basename(variant_set)
	
	mut <- read_parquet(variant_set)
	
	# calculate VAF and binomial test p-value for each mutation
	# NB  In calculating binom_p, it would be more accurate to consider
	#     copy-number states; e.g. in regions with CN of 3, the minor
	#     allele would only be 33% instead of 50%.
	mut <- mutate(mut,
		t_vaf = alt_ad / (alt_ad + ref_ad),
		binom_p = binom_tests(alt_ad, alt_ad + ref_ad, 0.5, alternative="less")
	)

	mut <- ungroup(mutate(
			group_by(mut, sample_name), 
			binom_q = p.adjust(binom_p, method = "BH")
		))

	# filter likely germline mutations from unmatched samples
	mut <- mutate(mut,
		pass = (binom_q < 1e-3 & alt_ad >= 5)
	)

	## Filter common variants based on population allele frequency
	## popaf is based on mutect2 annotations from gnomAD (gatk-best-practices)
	clin_thresh <- 0.0002 # 0.02%
	mut <- mutate(mut, common = (!is.na(popaf) & popaf >= clin_thresh))

	mut_clean <- mut[mut$pass & !mut$common, ]
	write_parquet(mut_clean, file.path(dataset, gsub(".parquet", "-gl.parquet", filename)))

	for (sample in unique(mut_clean$sample_name)){
		message(sprintf("	Writing filtered sample variants: %s", sample))
		mut_s <- mut_clean[mut_clean$sample_name == sample,  c("chrom", "pos", "ref", "alt")]

		outdir <- file.path("..", dataset, "filtered_pass-orientation-dp10-blacklist-gl", sample)
		dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
		write.table(mut_s, file.path(outdir, sprintf("%s.tsv", sample)), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
	}

}


