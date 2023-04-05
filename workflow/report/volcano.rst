Volcano plots visualizing differential expression analysis results of analysis {{snakemake.wildcards["analysis"]}}.
The feature's effect size as log2 fold change is indicated on the x-axis and significance as -log10 of {{snakemake.wildcards["pval_type"]}} on the y-axis. 
{{snakemake.wildcards["feature_list"]}} features are highlighted.