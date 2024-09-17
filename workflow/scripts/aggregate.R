#### load libraries & utility function 
library("ggplot2")
# library("patchwork")
library("data.table")
# library(reshape2)

# source utility functions
# source("workflow/scripts/utils.R")
snakemake@source("./utils.R")

# inputs
dea_result_path <- snakemake@input[["dea_results"]]

# outputs
dea_stats_path <- snakemake@output[["dea_stats"]]
all_features_path <- snakemake@output[["all_features"]]
dea_stats_plot_path <- snakemake@output[["dea_stats_plot"]]
dea_pvalue_plot_path <- snakemake@output[["dea_pvalue_plot"]]

# parameters
adj_pval <- as.numeric(snakemake@params[["adj_pval"]]) # 0.05
lfc <- as.numeric(snakemake@params[["lfc"]]) # 0
ave_expr <- as.numeric(snakemake@params[["ave_expr"]]) # 0
score_formula <- snakemake@params[["score_formula"]] # "-log10(dea_results$P.Value)*sign(dea_results$logFC)"

# plot specifications
width <- 0.25
height <- 5

# make result directory for feature lists if not exist
results_path <- dirname(file.path(all_features_path))
if (!dir.exists(results_path)){
    dir.create(results_path, recursive = TRUE)
}

### load DEA results
# dea_results <- read.csv(file=file.path(dea_result_path))
dea_results <- data.frame(fread(file.path(dea_result_path), header=TRUE))
groups <- unique(dea_results$group)

### save list of all expressed (and all significant features) for downstream analysis (eg as background in enrichment analyses)
all_features <- unique(dea_results$feature)
write(all_features, file.path(results_path, "ALL_features.txt"))

if("feature_name" %in% colnames(dea_results)){
    all_features_annot <- unique(dea_results$feature_name)
    all_features_annot <- all_features_annot[all_features_annot != ""]
    write(all_features_annot, file.path(results_path, "ALL_features_annot.txt"))
}

# write out as bed file if feature type is region
# (i.e. there is a column called "gencode_chr" in dea_results
# which is added by limma.R if "peak_id" is a col in consensus_annotations.csv)
bed_cols <- c("gencode_chr", "gencode_start", "gencode_end")
if ("gencode_chr" %in% colnames(dea_results)){
    all_features_bed <- unique(dea_results[, bed_cols])
    fwrite(all_features_bed, file=file.path(results_path, "ALL_features_annot.bed"), sep="\t", row.names=FALSE, col.names=FALSE)
}

# determine and save feature scores for each gene and group for downstream analysis e.g., preranked GSEA
if (score_formula!=""){
    dea_results$score <- eval(parse(text=score_formula))
    
    for (group in groups){
        tmp_features <- dea_results[(dea_results$group==group),c("feature","score")]
#         write.csv(tmp_features, file.path(results_path,paste0(group,"_featureScores.csv")), row.names=FALSE)
        fwrite(as.data.frame(tmp_features), file=file.path(results_path,paste0(group,"_featureScores.csv")), row.names=FALSE)


        if("feature_name" %in% colnames(dea_results)){
            tmp_features <- dea_results[(dea_results$group==group),c("feature_name","score")]
            tmp_features <- tmp_features[tmp_features$feature_name != "",]
#             write.csv(tmp_features, file.path(results_path,paste0(group,"_featureScores_annot.csv")), row.names=FALSE)
            fwrite(as.data.frame(tmp_features), file=file.path(results_path,paste0(group,"_featureScores_annot.csv")), row.names=FALSE)
        }
    }
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
                                             
# write.csv(dea_filtered_stats_df, file=file.path(dea_stats_path), row.names=TRUE)
fwrite(as.data.frame(dea_filtered_stats_df), file=file.path(dea_stats_path), row.names=TRUE)

### save differential feature lists from filtered DEA results for downstream analysis (eg enrichment analysis)
for (group in unique(dea_filtered_results$group)){
    for (direction in unique(dea_filtered_results$direction)){
        
        tmp_features <- dea_filtered_results[(dea_filtered_results$group==group) & (dea_filtered_results$direction==direction),"feature"]
        tmp_features <- unique(tmp_features)
        write(tmp_features, file.path(results_path, paste0(group,"_",direction,"_features.txt")))
    }
}

# write out as bed file if feature type is region
# (i.e. there is a column called "gencode_chr" in dea_results,
# which is added by limma.R if "peak_id" is a col in consensus_annotations.csv)
if("gencode_chr" %in% colnames(dea_results)){
    for (group in unique(dea_filtered_results$group)){
        for (direction in unique(dea_filtered_results$direction)){
            tmp_bed <- dea_filtered_results[(dea_filtered_results$group==group) & (dea_filtered_results$direction==direction), bed_cols]
            tmp_bed <- unique(tmp_bed)
            fwrite(tmp_bed, file=file.path(results_path, paste0(group,"_",direction,"_features.bed")), sep="\t", row.names=FALSE, col.names=FALSE)
        }
    }
}

# write out feature_names
if("feature_name" %in% colnames(dea_results)){
    for (group in unique(dea_filtered_results$group)){
        for (direction in unique(dea_filtered_results$direction)){
            tmp_features <- dea_filtered_results[(dea_filtered_results$group==group) & (dea_filtered_results$direction==direction),"feature_name"]
            tmp_features <- tmp_features[tmp_features != ""]
            tmp_features <- unique(tmp_features)
            write(tmp_features, file.path(results_path, paste0(group,"_",direction,"_features_annot.txt")))
        }
    }
}                                             
                                             
### visualize & save DEA statistics
width_panel <- length(groups) * width + 2
                                             
# dea_results_p <- ggplot(dea_filtered_results, aes(x=group, fill=direction)) + 
#                                              geom_bar() + 
#                                              xlab("groups") +
#                                              ylab("number of differential features") +
#                                              scale_y_continuous(limits = c(0, NA), expand = c(0, 0)) +
#                                              scale_x_discrete(label=addline_format) +
#                                              scale_fill_manual(values=list(up="red", down="blue"), drop=FALSE) +
#                                              custom_theme +
#                                              theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1, size = 6)) # rotates the x-Axis

dea_filtered_stats_df$total <- NULL
dea_filtered_stats_df[,"down"] <- -1 * dea_filtered_stats_df[,"down"]
plot_stats_df <- stack(dea_filtered_stats_df)
colnames(plot_stats_df) <- c("n_features","direction")
plot_stats_df$groups <- rep(rownames(dea_filtered_stats_df), ncol(dea_filtered_stats_df))

# plot
dea_filtered_results_p <- ggplot(plot_stats_df, aes(x=groups, y=n_features, fill=direction)) + 
                                             geom_bar(stat="identity", position="identity") +
                                             xlab("groups") +
                                             ylab("number of differential features") +
                                             scale_fill_manual(values=list("down"="blue", "up"="red"), drop=FALSE) +
                                             scale_y_continuous(labels = function(y) sapply(y, function(y) ifelse(y < 0, paste0(sub("-", "", as.character(y))), y))) +
                                             custom_theme +
                                             theme(#legend.position = "none",
                                                   axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 6)
                                                  )

                                             
# save plot
# options(repr.plot.width=width_panel, repr.plot.height=height)
# print(dea_results_p)
ggsave_new(filename = "stats", 
           results_path=dirname(dea_stats_plot_path), 
           plot=dea_filtered_results_p, #dea_results_p, 
           width=width_panel, 
           height=height)
                                             
                                             
# plot P-value distribution as sanity check
# pvalue_plots <- list()   
                                                                                            

                                             
for (group in unique(dea_results$group)){
    tmp_dea_results <- dea_results[dea_results$group==group, ]

#     pvalue_plots[[group]] <-
    pvalue_plot <- ggplot(tmp_dea_results, aes(x=P.Value, fill=factor(round(AveExpr)))) + geom_histogram(bins=100) + theme_bw(16) + ggtitle(addline_format(group)) + custom_theme
    
    ggsave_new(filename = group,
               results_path=dea_pvalue_plot_path, 
               plot=pvalue_plot, 
               width=4, 
               height=4
              )

}
                                             
# width <- 4
# height <- 4

# n_col <- min(10, length(unique(dea_results$group)))

# width_panel <- n_col * width + 1
# height_panel <- height * ceiling(length(unique(dea_results$group))/n_col)

# pvalue_plots_panel <- wrap_plots(pvalue_plots, ncol = n_col, guides = "collect")

# save plot
# options(repr.plot.width=width_panel, repr.plot.height=height_panel)
# print(pvalue_plots_panel)

# ggsave_new(filename = "pvalue_distribution", 
#            results_path=dirname(dea_stats_plot_path), 
#            plot=pvalue_plots_panel, 
#            width=width_panel, 
#            height=height_panel)