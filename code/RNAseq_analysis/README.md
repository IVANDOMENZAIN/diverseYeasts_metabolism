# Description
Scripts for performing gene differential expression analysis for datasets of transcripts reads for **C. intermedia** growing on several carbon sources (batch growth). The analysis assumes a negative binomial distribution of fold-changes for read counts between a given condition and reference samples (glucose). P-values are corrected for multiple testing using the Benjamini-Hochberg correction method.

# Instructions
Run the scripts in the following order.
 
1. `RNA_DE_analysis.m`: Perform DE analysis for growth on: celobiose, galactose, lactose and xylose.
2. `mapRNAdata2model.m`: Generate files including DE results and gene ids used in the model for facilitated mapping of results and data integration.

# Results
Results files can be found in the folder `results/RNA_DE_analysis`.
