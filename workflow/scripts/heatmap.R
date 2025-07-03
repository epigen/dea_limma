#### load libraries & utility function 
library("pheatmap")
library("ggplot2")
library("ggplotify")
library("reshape2")
library("data.table")

# source utility functions
# source("workflow/scripts/utils.R")
# snakemake@source("./utils.R") # does not work when loaded as module (https://github.com/snakemake/snakemake/issues/2205)
source(snakemake@params[["utils_path"]])

# inputs
dea_result_path <- snakemake@input[["dea_results"]]

# outputs
dea_lfc_heatmap_path <- snakemake@output[["dea_lfc_heatmap"]]

# parameters
adj_pval <- as.numeric(snakemake@params[["adj_pval"]]) # 0.05
lfc <- as.numeric(snakemake@params[["lfc"]]) # 0
ave_expr <- as.numeric(snakemake@params[["ave_expr"]]) # 0
feature_list_name <- snakemake@wildcards[["feature_list"]]

# plot specifications
width <- 0.15
height <- 0.15

### load DEA results
dea_results <- data.frame(fread(file.path(dea_result_path), header=TRUE))

# quit early and create empty result files, if there are no results
if(nrow(dea_results)==0){
    for (path in unlist(snakemake@output)) {
        dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE)
        file.create(path, showWarnings = FALSE)
    }
    quit(save = "no", status = 0)
}

# generate or load feature list
if (feature_list_name=="FILTERED"){
    feature_list <- unique(dea_results[(dea_results$adj.P.Val <= adj_pval) & 
                                  (abs(dea_results$logFC) >= lfc) & 
                                  (dea_results$AveExpr >= ave_expr),'feature'])
}else{
    feature_list_path <- snakemake@config[["feature_lists"]][[feature_list_name]]
    feature_list <- scan(file.path(feature_list_path), character())
    
    # if feature annotation is used then map annotation to features
    if("feature_name" %in% colnames(dea_results)){
        feature_list <- unique(dea_results[dea_results$feature_name %in% feature_list, 'feature'])
    }
}

# if no features in the results, end early with empty plot
if(length(feature_list)==0){
    ggsave_new(filename = feature_list_name, 
           results_path=dirname(dea_lfc_heatmap_path), 
           plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Features not found in DEA results.") + theme_void(), 
           width=4, 
           height=1)
    
    quit(save = "no", status = 0)
}

# subset DEA results according to feature list
dea_results <- dea_results[dea_results$feature %in% feature_list, ,drop = FALSE]

# make LFC dataframe
lfc_df <- reshape2::dcast(dea_results, feature ~ group, value.var = 'logFC')
rownames(lfc_df) <- lfc_df$feature
lfc_df$feature <- NULL

# make adjusted p-value dataframe
adjp_df <- reshape2::dcast(dea_results, feature ~ group, value.var = 'adj.P.Val')
rownames(adjp_df) <- adjp_df$feature
adjp_df$feature <- NULL

# set NA values to 0 (NA because below LFC threshold during testing or filtering -> can this still happen?)
lfc_df[is.na(lfc_df)] <- 0

# indicate significance
adjp_df[adjp_df<=adj_pval] <- "*"
adjp_df[adjp_df>adj_pval] <- ""
adjp_df[is.na(adjp_df)] <- "" # can this still happen?

### visualize LFC of DEA results as heatmap
height_panel <-  if (nrow(lfc_df)<100) (height * nrow(lfc_df) + 2) else 5
width_panel <- width * ncol(lfc_df) + 2

# format rownames for plotting
if("feature_name" %in% colnames(dea_results)){
    labels_row <- dea_results[match(rownames(lfc_df), dea_results$feature), 'feature_name']
}else{
    labels_row <- rownames(lfc_df)
}                                     

# make heatmap only if less than 50000 features
if(nrow(lfc_df)<50000){
    lfc_heatmap <- as.ggplot(pheatmap(lfc_df,
                                      display_numbers = if(nrow(lfc_df)<100) adjp_df else FALSE,
                                      main=paste0("logFC of ", feature_list_name," features"),
                                      cluster_cols = ifelse(ncol(lfc_df)>1, TRUE, FALSE),
                                      cluster_rows = ifelse(nrow(lfc_df)>1, TRUE, FALSE),
                                      show_rownames = ifelse(nrow(lfc_df)<100, TRUE, FALSE),
                                      labels_row = labels_row,
                                      show_colnames = TRUE,
                                      fontsize = 5,
                                      fontsize_number = 10,
                                      angle_col = 45,
                                      treeheight_row = 10,
                                      treeheight_col = 10,
                                      cellwidth = 10,
                                      cellheight = ifelse(nrow(lfc_df)<100, 10, NA),
                                      breaks=seq(-max(abs(lfc_df)), max(abs(lfc_df)), length.out=200),
                                      color=colorRampPalette(c("blue", "white", "red"))(200),
                                      annotation_names_col = F,
                                      silent = TRUE
                                     )
                            )
}else{
    lfc_heatmap <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Too many features to cluster and visualize.") + theme_void()
}

# save plot
# options(repr.plot.width=width_panel, repr.plot.height=height_panel)
# print(lfc_heatmap)

ggsave_new(filename = feature_list_name, 
           results_path=dirname(dea_lfc_heatmap_path), 
           plot=lfc_heatmap, 
           width=width_panel, 
           height=height_panel)
