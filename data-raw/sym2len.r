gene.lengths.fpath <- "/path/to/sym_gene_lengths.txt" # it has to be explained how it was created

read_gene_lengths <- function(gene.lengths.fpath) {
	gl <- read.delim(gene.lengths.fpath, sep = "\t", head = TRUE, comment.char="#")
	sym2length <- gl$Length
	names(sym2length) <- gl$Symbol
	sym2length
}

sym2len <- read_gene_lengths(gene.lengths.fpath)

setwd("/home/josephus/local/src-r")
usethis::create_project(path="rnaseqeasy")
usethis::use_data(sym2len, compress="gzip")
