#### load libraries & utility function 
library(ggplot2)

# source utility functions
source("workflow/scripts/utils.R")

# inputs
dea_result_path <- snakemake@input[["dea_results"]] #"/nobackup/lab_bock/projects/macroIC/results/CC001/dea_limma/CC001_thp1_filtered/DEA_results.csv" 

# outputs
dea_stats_path <- snakemake@output[["dea_stats"]] #"/nobackup/lab_bock/projects/macroIC/results/CC001/dea_limma/CC001_thp1_filtered/DEA_stats.csv"
dea_filtered_lfc_path <- snakemake@output[["dea_filtered_lfc"]] #"/nobackup/lab_bock/projects/macroIC/results/CC001/dea_limma/CC001_thp1_filtered/DEA_FILTERED_LFC.csv"
dea_stats_plot_path <- snakemake@output[["dea_stats_plot"]] #"/nobackup/lab_bock/projects/macroIC/results/CC001/dea_limma/CC001_thp1_filtered/plots/DEA_stats.png"


# parameters
adj_pval <- snakemake@params[["adj_pval"]] # 0.05
lfc <- snakemake@params[["lfc"]] # 0
ave_expr <- snakemake@params[["ave_expr"]] # 0

# plot specifications
width <- 0.5
height <- 5

# make result directory for feature lists if not exist
results_path <- file.path(dirname(file.path(dea_stats_path)),"feature_lists")
if (!dir.exists(results_path)){
    dir.create(results_path, recursive = TRUE)
}

### load DEA results
dea_results <- read.csv(file=file.path(dea_result_path))


### save list of all expressed features for downstream analysis (eg as background in enrichment analyses)
all_features <- unique(dea_results$genes)
write(all_features, file.path(results_path, "ALL_features.txt"))


# annotate differential direction (up or down)
dea_results$direction <- as.character(lapply(dea_results$logFC, function(x) if(x>0){"up"}else{"down"}))

### aggregate & save FILTERED DEA statistics
dea_filtered_results <- dea_results[(dea_results$adj.P.Val <= adj_pval) & 
                                  (abs(dea_results$logFC) >= lfc) & 
                                  (dea_results$AveExpr >= ave_expr), ]
                                                                                          
dea_filtered_stats <- table(dea_filtered_results$group, dea_filtered_results$direction)
dea_filtered_stats_df <- as.data.frame.matrix(dea_filtered_stats)
dea_filtered_stats_df$total <- rowSums(dea_filtered_stats_df)
                                             
write.csv(dea_filtered_stats_df, file=file.path(dea_stats_path), row.names=TRUE)

### aggregate & save LFC matrix from filtered DEA result features
groups <- unique(dea_results$group)
features <- unique(dea_filtered_results$genes)
                                             
lfc_df <- data.frame(matrix(nrow=length(features), ncol=length(groups), dimnames=list(features, groups)))

for (group in unique(dea_results$group)){
    tmp_dea_results <- dea_results[dea_results$group==group, ]
    rownames(tmp_dea_results) <- tmp_dea_results$genes
    lfc_df[features, group] <- tmp_dea_results[features, 'logFC']
}

write.csv(lfc_df, file=file.path(dea_filtered_lfc_path), row.names=TRUE)

### save differential feature lists from filtered DEA results for downstream analysis                                             
for (group in unique(dea_filtered_results$group)){
    for (direction in unique(dea_filtered_results$direction)){
        
        tmp_features <- dea_filtered_results[(dea_filtered_results$group==group) & (dea_filtered_results$direction==direction),"genes"]
        write(tmp_features, file.path(results_path,paste0(group,"_",direction,"_features.txt")))
    }
}
                                             
### visualize & save DEA statistics
width_panel <- length(groups) * width + 1
                                             
dea_results_p <- ggplot(dea_filtered_results, aes(x=group, fill=direction)) + 
                                             geom_bar() + 
                                             xlab("groups") +
                                             ylab("number of differential features") +
                                             scale_y_continuous(limits = c(0, NA), expand = c(0, 0)) +
                                             scale_fill_manual(values=list(up="red", down="blue"), drop=FALSE) +
                                             custom_theme +
                                             theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1)) # rotates the x-Axis
                                             
                                             
# save plot
# options(repr.plot.width=width_panel, repr.plot.height=height)
# print(dea_results_p)
ggsave_new(filename = "DEA_stats", 
           results_path=dirname(dea_stats_plot_path), 
           plot=dea_results_p, 
           width=width_panel, 
           height=height)