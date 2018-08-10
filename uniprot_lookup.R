#!/usr/bin/env Rscript
library(UniProt.ws)
library(clusterProfiler)
library(argparse)

setwd("/home/matthew/lab_root/pathway_jamb")

# Single lookup (meant for collecting results for amigious gene IDs e.g. multiple IDs)
lookup_gene <- function(gene_id) {
		results <- select(up, keys=gene_id, columns=c("UNIPROTKB", "DATABASE(PFAM)", "GO-ID"), keytype="ENSEMBL_GENOMES")
		print(results)
}

# Often there are multiple uniprot entries for a given gene name
# This cleans out some of the lower quality entries if possible
clean_results <- function(results) {
		# if there is only one result anyway, there is no need to attempt to clean anything
		if (dim(results)[[1]] <= 1) {
			return(results)
		}
	
		## There are some entries that are completely NA except for gene name & uniprot id
		## They seem to correspond to indica entries, let's remove these rows
		results <- results[!(is.na(results$GO) & is.na(results$PFAM)),]
		
		# If we have multiple uniprot IDs one greater than 6 characters, one shorter
		# then pick the 6 digit one only
		# TODO: It seems like we are missing some good entries when using the API that
		# show up when using the web service at uniprot.org
		ids_size <- nchar(results$UNIPROTKB)
		if (6 %in% ids_size) {
			results <- results[ids_size == 6,]
		}
	
		# Pick out the highly rated curated entries	
		# if score is 5 out of 5, pick that one
		if ("5 out of 5" %in% results$SCORE) {
			results <- results[results$SCORE == "5 out of 5",]
		}
	return(results)
}

get_uniprot_ids <- function(gene_ids, up) {
	uniprot_ids <- list()
	GOs <- list()
	pfams <- list()
	for (i in 1:length(gene_ids)) {
		# don't even try to lookup an empty gene_id
		if (gene_ids[[i]] == "") {
			uniprot_ids[[i]] <- ""
			GOs[[i]] <- ""
			pfams[[i]] <- ""
			next
		}
		print(paste("looking up:", gene_ids[[i]]))
		# TODO: can we query all gene_ids at once while still checking each for errors?
		results <- tryCatch(
			{
				results <- select(up, keys=gene_ids[[i]], columns=c("UNIPROTKB", "DATABASE(PFAM)", "GO-ID", "SCORE"), keytype="ENSEMBL_GENOMES")
			},
			error = function(cond) {
				print(paste("Error accessing Uniprot database: ", cond))
				return(NULL)
			}
		)
		if (is.null(results)) {
			uniprot_ids[[i]] <- ""
			GOs[[i]] <- ""
			pfams[[i]] <- ""
			next;
		}
		
		names(results) <- c("GENE", "UNIPROTKB", "PFAM", "GO", "SCORE")
		results <- clean_results(results)
		print(results)
		
		if (dim(results)[[1]] > 1) {
			print(paste("Got more than one result for", gene_ids[[i]], "needs manual curation."))
			uniprot_ids[[i]] <- ""
			GOs[[i]] <- ""
			pfams[[i]] <- ""
		} else if (dim(results)[[1]] == 1) {
			uniprot_ids[[i]] <- results$UNIPROTKB	
			GOs[[i]] <- results$GO
			pfams[[i]] <- results$PFAM
		} else {
			print(paste("Couldn't find", gene_ids[[i]], "in UniProt."))
			uniprot_ids[[i]] <- ""
			GOs[[i]] <- ""
			pfams[[i]] <- ""
		}
	}
	
	# TODO: Investigate/remove terminal ';' in PFAMs
	return (list(uniprot_ids=uniprot_ids, GOs=GOs, pfams=pfams))
}
main <- function(input_file, output_file, taxId) {
	up <- UniProt.ws::UniProt.ws(taxId=taxId)
	input_table <- read.csv(input_file, header=TRUE, sep="\t")
	#names(input_table) <- c("GENE_NAMES", "MSU_ID", "RAP_ID", "UNIPROT_ID", "PFAM", "SUBCELLULAR", "TRANSMEMBRANE_DOMAIN","GO", "EC#", "References (PMID)", "Curator Name",	"Remark by Curator", "Associated Reaction", "Associated pathway")
	gene_ids <- as.character(input_table$RAP_geneID)
	uniprot_results <- get_uniprot_ids(gene_ids, up)
	input_table$Uniprot_ID <- unlist(uniprot_results$uniprot_ids)
	input_table$Pfam_domain <- unlist(uniprot_results$pfams)
	input_table$GO <- unlist(uniprot_results$GOs)
	write.table(input_table, file=output_file, quote=FALSE, sep='\t', row.names=FALSE)
}

parser <- ArgumentParser()
parser$add_argument("-i", "--input", type="character", dest="input", help="Input file to parse(TSV format)")
parser$add_argument("-o", "--output", type="character", dest="output", help="Output file name")
parser$add_argument("--taxId", type="integer", dest="taxID", default=39947, help="taxId(default 39947)")
args <- parser$parse_args()
if (is.null(args$input)) {
	cat("Missing required input file arugment(specify with -i or --input)\n")
} else if (is.null(args$output)) {
	cat("Missing output file argument(specify with -o or --output)\n")
} else {
	main(args$input, args$output, args$taxID)
}