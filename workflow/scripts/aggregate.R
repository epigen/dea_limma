#### load libraries & utility function 
library(ggplot2)
library(patchwork)
library(ggtext)

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
all_features <- unique(dea_results$feature)
write(all_features, file.path(results_path, "ALL_features.txt"))

if("feature_name" %in% colnames(dea_results)){
    all_features_annot <- unique(dea_results$feature_name)
    all_features_annot <- all_features_annot[all_features_annot != ""]
    write(all_features_annot, file.path(results_path, "ALL_features_annot.txt"))
}

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
features <- unique(dea_filtered_results$feature)
                                             
lfc_df <- data.frame(matrix(nrow=length(features), ncol=length(groups), dimnames=list(features, groups)))

for (group in unique(dea_results$group)){
    tmp_dea_results <- dea_results[dea_results$group==group, ]
    rownames(tmp_dea_results) <- tmp_dea_results$feature
    lfc_df[features, group] <- tmp_dea_results[features, 'logFC']
}

write.csv(lfc_df, file=file.path(dea_filtered_lfc_path), row.names=TRUE)
                                             
if("feature_name" %in% colnames(dea_results)){
    lfc_df$feature_name <- tmp_dea_results[rownames(lfc_df), 'feature_name']
    lfc_df <- lfc_df[,c(ncol(lfc_df),1:(ncol(lfc_df)-1))]
    write.csv(lfc_df, file=file.path(dirname(file.path(dea_stats_path)),"DEA_FILTERED_LFC_annot.csv"), row.names=FALSE)
}

### save differential feature lists from filtered DEA results for downstream analysis (eg enrichment analysis)
for (group in unique(dea_filtered_results$group)){
    for (direction in unique(dea_filtered_results$direction)){
        
        tmp_features <- dea_filtered_results[(dea_filtered_results$group==group) & (dea_filtered_results$direction==direction),"feature"]
        write(tmp_features, file.path(results_path,paste0(group,"_",direction,"_features.txt")))
    }
}
                                             
if("feature_name" %in% colnames(dea_results)){
    for (group in unique(dea_filtered_results$group)){
        for (direction in unique(dea_filtered_results$direction)){

            tmp_features <- dea_filtered_results[(dea_filtered_results$group==group) & (dea_filtered_results$direction==direction),"feature_name"]
            tmp_features <- tmp_features[tmp_features != ""]
            write(tmp_features, file.path(results_path,paste0(group,"_",direction,"_features_annot.txt")))
        }
    }
}                                             
                                             
### visualize & save DEA statistics
width_panel <- length(groups) * width + 1
                                             
dea_results_p <- ggplot(dea_filtered_results, aes(x=group, fill=direction)) + 
                                             geom_bar() + 
                                             xlab("groups") +
                                             ylab("number of differential features") +
                                             scale_y_continuous(limits = c(0, NA), expand = c(0, 0)) +
                                             scale_x_discrete(label=addline_format) +
                                             scale_fill_manual(values=list(up="red", down="blue"), drop=FALSE) +
                                             custom_theme +
                                             theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1, size = 6)) # rotates the x-Axis
                                             
                                             
# save plot
# options(repr.plot.width=width_panel, repr.plot.height=height)
# print(dea_results_p)
ggsave_new(filename = "DEA_stats", 
           results_path=dirname(dea_stats_plot_path), 
           plot=dea_results_p, 
           width=width_panel, 
           height=height)
                                             
                                             
# plot P-value distribution as sanity check
pvalue_plots <- list()                                             
                                             
for (group in unique(dea_results$group)){
    tmp_dea_results <- dea_results[dea_results$group==group, ]

    pvalue_plots[[group]] <- ggplot(tmp_dea_results, aes(x=P.Value, fill=factor(round(AveExpr)))) + geom_histogram(bins=100) + theme_bw(16) + ggtitle(addline_format(group)) + custom_theme
}
                                             
width <- 4
height <- 4

n_col <- min(10, length(unique(dea_results$group)))

width_panel <- n_col * width + 1
height_panel <- height * ceiling(length(unique(dea_results$group))/n_col)

pvalue_plots_panel <- wrap_plots(pvalue_plots, ncol = n_col, guides = "collect")

# save plot
# options(repr.plot.width=width_panel, repr.plot.height=height_panel)
# print(pvalue_plots_panel)

ggsave_new(filename = "DEA_pvalue_distribution", 
           results_path=dirname(dea_stats_plot_path), 
           plot=pvalue_plots_panel, 
           width=width_panel, 
           height=height_panel)