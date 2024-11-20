#### load libraries & utility function 
library("ggplot2")
library("data.table")

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

# determine and save feature scores for each gene and group for downstream analysis e.g., preranked GSEA
if (score_formula!=""){
    dea_results$score <- eval(parse(text=score_formula))
    
    for (group in groups){
        tmp_features <- dea_results[(dea_results$group==group),c("feature","score")]
        fwrite(as.data.frame(tmp_features), file=file.path(results_path,paste0(group,"_featureScores.csv")), row.names=FALSE)


        if("feature_name" %in% colnames(dea_results)){
            tmp_features <- dea_results[(dea_results$group==group),c("feature_name","score")]
            tmp_features <- tmp_features[tmp_features$feature_name != "",]
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

# Handling the exception when no significant genes are found in either up or down categories
if (!all(c("up", "down") %in% colnames(dea_filtered_stats_df))) {
    dea_filtered_stats_df[setdiff(c("up", "down"), colnames(dea_filtered_stats_df))] <- 0
}

dea_filtered_stats_df$total <- rowSums(dea_filtered_stats_df)
fwrite(as.data.frame(dea_filtered_stats_df), file=file.path(dea_stats_path), row.names=TRUE)

### save differential feature lists from filtered DEA results for downstream analysis (eg enrichment analysis)
for (group in groups){
    for (direction in c("up", "down")){
# for (group in unique(dea_filtered_results$group)){
    # for (direction in unique(dea_filtered_results$direction)){
        tmp_features <- dea_filtered_results[(dea_filtered_results$group==group) & (dea_filtered_results$direction==direction),"feature"]
        tmp_features <- unique(tmp_features)
        write(tmp_features, file.path(results_path,paste0(group,"_",direction,"_features.txt")))

        if("feature_name" %in% colnames(dea_results)){
            tmp_features <- dea_filtered_results[(dea_filtered_results$group==group) & (dea_filtered_results$direction==direction),"feature_name"]
            tmp_features <- tmp_features[tmp_features != ""]
            tmp_features <- unique(tmp_features)
            write(tmp_features, file.path(results_path,paste0(group,"_",direction,"_features_annot.txt")))
        }
    }
}
                                             
if (nrow(dea_filtered_stats) != 0) { # If there are DEA results after filtering
    ### visualize DEA statistics
    groups <- unique(dea_filtered_results$group)
    width_panel <- length(groups) * width + 2

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
                                             # geom_text(aes(label=abs(n_features), vjust=ifelse(n_features < 0, 1.5, -0.5), hjust=0.5), size=2) +
                                             geom_text(aes(label=ifelse(n_features == 0, '', abs(n_features)), vjust=ifelse(n_features < 0, 1.5, -0.5), hjust=0.5), size=2) +
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
                                             
} else { # If there are no filtered features.
	# empty stats plot
	ggsave_new(filename = "stats", 
           results_path=dirname(dea_stats_plot_path), 
           plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No filtered features in DEA results.") + theme_void(), 
           width=4, 
           height=1)
}
                                             
# plot P-value distribution as sanity check                                             
for (group in unique(dea_results$group)){
    tmp_dea_results <- dea_results[dea_results$group==group, ]

    pvalue_plot <- ggplot(tmp_dea_results, aes(x=P.Value, fill=factor(round(AveExpr)))) + 
        geom_histogram(bins=100) + 
        theme_bw(16) + 
        ggtitle(addline_format(group)) + 
        custom_theme
    
    ggsave_new(filename = group,
               results_path=dea_pvalue_plot_path, 
               plot=pvalue_plot, 
               width=4, 
               height=4
              )

}
