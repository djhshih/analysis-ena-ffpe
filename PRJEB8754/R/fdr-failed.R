#!/usr/bin/env Rscript
# Identify entries that fail adaptive FDR threshold

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


library(argparser)
library(stringr)

p <- arg_parser("Select variants based on adaptive FDR")

p <- add_argument(p, "snvf", help="input SNV file")
p <- add_argument(p, "--out", help="output file name")
p <- add_argument(p, "--fp-cut", help="expected false positive threshold", default=0.5)

argv <- parse_args(p)

if (is.na(argv$out)) {
	out.fn <- paste(str_extract(argv$snvf, ".+_"), "failed.snv", sep="")
} else {
	out.fn <- argv$out;
}

# The input may be an empty file:
fileinfo <- file.info(argv$snvf)

if(fileinfo$size == 0 ){
	write.table(data.frame(), file=argv$out, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
}else{
	x <- read.table(argv$snvf, sep="\t", fill=TRUE, header=TRUE)

	# delete the NAs
	x <- x[complete.cases(x$FOBP), ]

	# set the 0 to machine precision
	x$FOBP[x$FOBP == 0] <- .Machine$double.eps

	# adjust p value
	x$q <- p.adjust(x$FOBP, "BH")

	# Use adaptive FDR to identify called SNVs
	pass <- adaptive_fdr_cut(x$q, argv$fp_cut)

	# identify all the failed SNVs in the SNV file
	x.failed <- x[!pass, ]

	write.table(x.failed, file=argv$out, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

}


