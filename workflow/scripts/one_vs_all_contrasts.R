#### load libraries
library("limma")
library("data.table")

#### configs
# input
lmfit_object_path <- snakemake@input[["lmfit_object"]]
model_matrix_path <- snakemake@input[["model_matrix"]]
feature_annotation_path <- snakemake@input[["feature_annotation"]]
metadata_path <- snakemake@input[["metadata"]]

# output
contrast_result_path <- snakemake@output[["contrast_results"]]
contrast_object_path <- snakemake@output[["contrast_object"]]
contrast_matrix_path <- snakemake@output[["contrast_matrix"]]

# parameters
ova_var <- snakemake@wildcards[["ova_var"]]
feature_annotation_col <- base::make.names(snakemake@params[["feature_annotation_col"]])[1]
eBayes_flag <- snakemake@params[["eBayes"]] #1
limma_trend <- snakemake@params[["limma_trend"]] #0

#### load data (model, design and feature_annotation)
fit <- readRDS(lmfit_object_path)
design <- data.frame(fread(file.path(model_matrix_path), header=TRUE), row.names=1)

### load feature annotation file (optional)
if (feature_annotation_path!=""){
    feature_annot <- data.frame(fread(file.path(feature_annotation_path), header=TRUE), row.names=1)
    print("feature_annot")
    print(dim(feature_annot))
}

#### identify groups for one-vs-all contrasts using the metadata table, column prefix & design matrix
meta <- data.frame(fread(file.path(metadata_path), header=TRUE), row.names=1)
if (grepl(':', ova_var)){
    # for interaction terms, need to generate all combinations of the levels
    # there could be multiple colons in the variable name for high-order interaction terms
    ova_vars <- unlist(strsplit(ova_var, ':'))
    individual_group_levels <- lapply(ova_vars, function(var) {
        unique(meta[, var])
    })
    # get combinations of all group levels
    level_combinations <- expand.grid(individual_group_levels, stringsAsFactors = FALSE)
    group_names <- apply(level_combinations, 1, function(row) {
        paste(row, collapse = ":")
    })
    ## add prefix to match the design matrix
    group_cols <- apply(level_combinations, 1, function(row) {
        paste0(ova_vars, row, collapse = ":")
    })
} else {
    group_names <- unique(meta[,ova_var])
    # add prefix to match the design matrix
    group_cols <- paste0(ova_var, group_names)
}

#### define contrasts
contrasts_all <- list()

for (i in seq(length(group_names))){
    group_name <- group_names[[i]]
    group_col <- group_cols[[i]]
    
    # list of all group_cols without current group group_col
    group_cols_wo_gr <- group_cols[group_cols != group_col]
    # define contrast for one-vs-all analysis
    contrasts_all[[group_name]] <- paste0(group_col, " - (", paste(group_cols_wo_gr, collapse=" + "), ") / ", length(group_cols_wo_gr))
    print(contrasts_all[[group_name]])
}

# generate contrast matrix based on contrast formulas and design matrix
contrast_matrix <- makeContrasts(contrasts=contrasts_all, levels = design)
# save contrast matrix
fwrite(as.data.frame(contrast_matrix), file=file.path(contrast_matrix_path), row.names=TRUE)

#### fit contrasts
fit2 <- contrasts.fit(fit, contrast_matrix)
# save fitted model object for manual downstream analyses
saveRDS(fit2, file = file.path(contrast_object_path))

# estimate/correct variance with eBayes (optional) without or with limma-trend (optional)
# fit2 <- eBayes(fit2)
if (eBayes_flag!=0){
    if(limma_trend==0){
        fit2 <- eBayes(fit2, robust=TRUE, trend=FALSE)
    }else{
        fit2 <- eBayes(fit2, robust=TRUE, trend=TRUE)
    }
}else{ # Warning: not recommended by limma author Gordon Smyth (https://support.bioconductor.org/p/35174/)
    # determine ordinary t statistic, p-value and B-statistic (lods) according to
    # limma userguide chapter 13.2 Fitted Model Objects: http://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf
    # bioconductor forum: https://support.bioconductor.org/p/113833/
    # bioconductor forum: https://support.bioconductor.org/p/35159/
    # note: B-statistic determined using eBayes, assuming it will not be used downstream
    
    fit2$t <- fit2$coef/fit2$stdev.unscaled/fit2$sigma
    fit2$p.value <- 2 * pt(-abs(fit2$t), df = fit2$df.residual)
    fit2$lods <- eBayes(fit2, robust=TRUE, trend=FALSE)$lods
}

# extract results
contrast_result <- data.frame()

for(coefx in colnames(coef(fit2))){
    tmp_res <- topTable(fit2, coef=coefx, number=nrow(fit2), sort.by="P")
    tmp_res$feature <- rownames(tmp_res)
    tmp_res <- tmp_res[,c(ncol(tmp_res),1:(ncol(tmp_res)-1))]
    rownames(tmp_res) <- NULL
    
    # beautify group names by replacing them (the contrast formula) with the names of the comparison
    tmp_res$group <- names(contrasts_all)[contrasts_all == coefx]

    # add feature_name
    if (feature_annotation_path!=""){
        tmp_res$feature_name <- feature_annot[tmp_res$feature, feature_annotation_col]
    }
    
    if(dim(contrast_result)[1]==0){
        contrast_result <- tmp_res
    }else{
        contrast_result <- rbind(contrast_result, tmp_res)
    }
}

# remove rows with adj.P.Val=NA
contrast_result <- contrast_result[!is.na(contrast_result$adj.P.Val),]

#### save results
fwrite(as.data.frame(contrast_result), file=file.path(contrast_result_path), row.names=FALSE)
