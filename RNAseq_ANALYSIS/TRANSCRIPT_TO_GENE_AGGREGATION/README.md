RNA-seq Transcript to Gene Count Aggregation
This R script processes RNA-seq transcript-level count data produced with salmon and aggregates it to gene-level counts. 

Requirements
tximport, readr, dplyr – all available in R.

Input Files
-	tx2gene.txt – produced via Biomart (a 2 column mapping file: transcript_ID, gene_ID).
-	counts_matrix.sf – produced via salmon as in (Saunderson and Dobson, 2025) (tab separated with sample columns).

What it does
1.	Loads transcript counts from counts_matrix.sf.
2.	Maps transcripts to genes using the tx2gene.txt mapping file.
3.	Aggregates transcript counts by summing all transcripts belonging to the same gene.
4.	Filters out unmapped transcripts.
5.	Exports gene-level counts to output

Output
“gene_level_counts.txt” - Tab-separated gene-level count matrix

Usage
Simply replace the example working directory with your own location which contains the transcript level “counts_matrix.sf” output from salmon, and run the script!


