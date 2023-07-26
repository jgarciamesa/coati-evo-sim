# Get sequences from Ensembl

library(httr)
library(jsonlite)
library(xml2)
library(seqinr)
library(stringr)

get_canonical_transcript = function(request) {
	r = fromJSON(toJSON(content(request)))
	return(r$Transcript[["id"]][r$Transcript$is_canonical==1])
}

getseq_main = function(gene_table, gene_id, out_folder) {
	server = "https://rest.ensembl.org/"

	genes = read.table(gene_table, header = FALSE, sep = "\t",
	                   stringsAsFactors=FALSE)
	# find row with gene_id
	download_ids = genes[which(genes == gene_id),]

	sequences = c()

	for(id in download_ids) {
		req = GET(paste(server, "lookup/id/", id, "?expand=1", sep = ""))
		stop_for_status(req)

		id = get_canonical_transcript(req)

		req = GET(paste(server, "sequence/id/", id, "?type=cds", sep = ""),
			content_type("text/x-fasta"))
		stop_for_status(req)

		sequences = append(sequences, content(req))
	}
	write.table(sequences,
	            file = paste(out_folder, "/", gene_id, ".fasta", sep=""),
	            quote = FALSE, row.names = FALSE, col.names = FALSE)
}

if(!interactive()) {
	ARGS = commandArgs(trailing=TRUE)
	getseq_main(gene_table = ARGS[1], gene_id = ARGS[2], out_folder = ARGS[3])
}
