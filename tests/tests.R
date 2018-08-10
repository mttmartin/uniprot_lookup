#!/usr/bin/env Rscript
source("../uniprot_lookup.R")

context("Testing result filtering")

up <- UniProt.ws::UniProt.ws(taxId=39947)

test_that('Mostly null entries removed', {
		results <- data.frame(c("Os06g0644200", "Os06g0644200"), c("Q67WN4", "Q67WN4"), c(NA, "PF00249"), c(NA, "GO:0003677"), c("1 out of 5","1 out of 5"), stringsAsFactors=FALSE)
		names(results) <- c("GENE", "UNIPROTKB", "PFAM", "GO", "SCORE")
		cleaned_results <- clean_results(results)
		expect_equal(nrow(cleaned_results), 1)
})

test_that('Six digit Uniprot IDs preferred', {
	results <- data.frame(c("Os06g0644200", "Os06g0644200"), c("Q67WN4", "A0A0P0XHA2"), c("PF00249", "PF00249"), c("GO:0003677", "GO:0003677"), c("1 out of 5","1 out of 5"), stringsAsFactors=FALSE)	
	names(results) <- c("GENE", "UNIPROTKB", "PFAM", "GO", "SCORE")
	cleaned_results <- clean_results(results)
	expect_equal(nrow(cleaned_results), 1)
})


test_that('Entries with highest scores are preferred', {
	results <- data.frame(c("Os06g0644200", "Os06g0644200"), c("Q67WN4","Q67WN4"), c("PF00249", "PF00249"), c("GO:0003677", "GO:0003677"), c("5 out of 5","1 out of 5"), stringsAsFactors=FALSE)	
	names(results) <- c("GENE", "UNIPROTKB", "PFAM", "GO", "SCORE")
	cleaned_results <- clean_results(results)
	expect_equal(nrow(cleaned_results), 1)
})

context("Testing Uniprot database lookup")

test_that('Empty gene ID works', {
		gene_id <- c("")
		results <- get_uniprot_ids(gene_id, up)
		expect_equal(length(results$uniprot_ids), 1)
		expect_equal(results$uniprot_ids[[1]], "")
})

test_that('Invalid gene ID works', {
		gene_id <- c("foobar")
		results <- get_uniprot_ids(gene_id, up)
		expect_equal(length(results$uniprot_ids), 1)
		expect_equal(results$uniprot_ids[[1]], "")
})

test_that('Valid gene ID works', {
		gene_id <- c("Os06g0614100")
		results <- get_uniprot_ids(gene_id, up)
		expect_equal(length(results$uniprot_ids), 1)
		#expect neq id ""
})

