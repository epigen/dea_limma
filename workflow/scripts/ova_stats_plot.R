#### load libraries & utility function 
library("ggplot2")
library("data.table")

# source utility functions
# source("workflow/scripts/utils.R")
# snakemake@source("./utils.R") # does not work when loaded as module (https://github.com/snakemake/snakemake/issues/2205)
source(snakemake@params[["utils_path"]])

# inputs
dea_result_paths <- snakemake@input[["dea_results"]]

# outputs
dea_stats_path <- snakemake@output[["dea_stats_OvA"]]
dea_stats_plot_path <- snakemake@output[["dea_stats_plot_OvA"]]

# parameters
adj_pval <- as.numeric(snakemake@params[["adj_pval"]]) # 0.05
lfc <- as.numeric(snakemake@params[["lfc"]]) # 0
ave_expr <- as.numeric(snakemake@params[["ave_expr"]]) # 0
score_formula <- snakemake@params[["score_formula"]] # "-log10(dea_results$P.Value)*sign(dea_results$logFC)"

# plot specifications
width <- 0.25
height <- 5

### load DEA results
dea_results <- data.frame()
for (path in dea_result_paths) {
    new_dea_results <- data.frame(fread(file.path(path), header=TRUE))
    # results for OvA are saved without group prefix, add this here
    group_name <- gsub(".*dea_limma/.*_OvA_([^/]+)/results.csv", "\\1", path)
    if (!any(grepl(group_name, new_dea_results$group))) {
        new_dea_results$group <- paste0(group_name, "_", new_dea_results$group)
    }
    dea_results <- rbind(dea_results, new_dea_results)
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
    missing_cols <- setdiff(c("up", "down"), colnames(dea_filtered_stats_df))
    if (length(missing_cols) == 2) {
        # just make a new dataframe, otherwise index fails completely
        dea_filtered_stats_df <- data.frame(up=0, down=0)
    } else {
        dea_filtered_stats_df[missing_cols] <- 0
    }
}

dea_filtered_stats_df$total <- rowSums(dea_filtered_stats_df)
fwrite(as.data.frame(dea_filtered_stats_df), file=file.path(dea_stats_path), row.names=TRUE)

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
                                             scale_fill_manual(values=list("up"="red", "down"="blue"), drop=FALSE) +
                                             scale_y_continuous(labels = function(y) sapply(y, function(y) ifelse(y < 0, paste0(sub("-", "", as.character(y))), y))) +
                                             geom_text(aes(label=ifelse(n_features == 0, '', abs(n_features)), vjust=ifelse(n_features < 0, 1.5, -0.5), hjust=0.5), size=2) +
                                             custom_theme +
                                             theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 6))

    # save plot
    # options(repr.plot.width=width_panel, repr.plot.height=height)
    # print(dea_results_p)
    ggsave_new(filename = "stats_OvA", 
           results_path=dirname(dea_stats_plot_path), 
           plot=dea_filtered_results_p, #dea_results_p, 
           width=width_panel, 
           height=height)
                                             
} else { # If there are no filtered features.
	# empty stats plot
	ggsave_new(filename = "stats_OvA", 
           results_path=dirname(dea_stats_plot_path), 
           plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No filtered features in DEA results.") + theme_void(), 
           width=4, 
           height=1)
}