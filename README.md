# Description
This program might be useful for automating some aspects of manual pathway curation. It was developed as part of the Pathway Curation Workshop at the International Conference on Biological Ontology as a prototype to examine automation of some preliminary steps for pathway curation. 

During the initial research steps for curation, a number of manual database lookups are required for each gene of interest. This program automates these lookup for simple cases not requiring human intervention. For more complicated scenarios, it simply leaves the entries blank for humans to manually handle.

# Dependencies
The following libraries are required
* Uniprot.ws
* clusterProfiler
* argparse

# Usage

The input file is expected to contain the following columns(with these exact column names): 

	* RAP_geneID
	* Pfam_domain
	* Uniprot_ID
	* GO

The program can be executed using:

```bash
./test_uniprot_lookup.R -i input_file.tsv -o output_file.tsv
```


By default the program looks for entries from Oryza sativa subsp. japonica, but a different NCBI taxonomic identifier can be specified using the experimental --taxId parameter.
