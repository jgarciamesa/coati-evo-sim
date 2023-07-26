# Get gene ID list from Ensembl

library(biomaRt)
library(dplyr, warn.conflicts = FALSE)

# As explained in bioconductor.org, if R version is 3.6 or above:
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
#
# BiocManager::install("biomaRt")
#
# For documentation of this version:
# browseVignettes("biomaRT")

get_geneId = function(output) {
	dir.create(file.path("data"), showWarnings = FALSE,
	           recursive = TRUE)

	# Query to Ensembl and select homo sapiens gene dataset
	ensembl_h = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
	# Filters to use
	filters = c("with_ccds")
	# Attributes to retrieve
	attributes = c("ensembl_gene_id")
	values = list(TRUE)#, TRUE)

	genes = getBM(attributes, filters, values, mart = ensembl_h, uniqueRows = TRUE)

	# sort
	genes = genes[order(genes$ensembl_gene_id),]

	write.table(genes, file = output, quote = FALSE, sep = "\t",
	            row.names = FALSE, col.names = FALSE)
}

if(!interactive()) {
	ARGS = commandArgs(trailingOnly = TRUE)
	get_geneId(output = ARGS[1])
}

