#!/usr/bin/env Rscript
library(argparser)
library(stringr)
library(glue)


# Identify entries that fail adaptive FDR threshold
# fp-cut depends on sequencing depth
# from https://github.com/djhshih/analysis-tcga-ffpe/blob/master/gdc/vcf-ffpe-snvf/analyze.R

adaptive_fdr_cut <- function(q, fp.cut) {
	n <- length(q);
	# adaptive threshold chosen to expect less than one false positive
	idx <- order(q);
	under <- which((1:n) * q[idx] < fp.cut);
	if (length(under) > 0) {
		top <- under[length(under)];
		pred <- logical(n);
		# select the top significant results
		pred[idx[1:top]] <- TRUE;
		pred
	} else {
		# none passes
		rep(FALSE, n)
	}
}



paths <- list.files("../../ffpe-snvf/mobsnvf/", pattern = ".mobsnvf.snv$", recursive = TRUE, full.names = TRUE)


for (path in paths) {

	x <- read.table(path, sep="\t", fill=TRUE, header=TRUE)

	y <- x[!complete.cases(x$FOBP), ]

	# delete the NAs
	x <- x[complete.cases(x$FOBP), ]

	# set the 0 to machine precision
	x$FOBP[x$FOBP == 0] <- .Machine$double.eps

	# adjust p value
	x$q <- p.adjust(x$FOBP, "BH")

	# Use adaptive FDR to identify called SNVs
	pass <- adaptive_fdr_cut(x$q, 0.5)

	# identify all the failed SNVs in the SNV file
	x.failed <- x[!pass, ]
	x.passed <- x[pass, ]

	if (dim(y)[1] != 0){
		y$q <- NA
		x.passed <- rbind(x.passed, y)
	}

	x.passed <- x.passed[order(x.passed$chrom, x.passed$pos, x.passed$ref, x.passed$alt), ]

	out_prefix <- file.path(dirname(path), tools::file_path_sans_ext(basename(path)))

	write.table(x.failed, file=paste0(out_prefix, ".artifacts.snv"), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
	write.table(x.passed, file=paste0(out_prefix, ".filtered.snv"), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

}




