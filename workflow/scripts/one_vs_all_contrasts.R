#### load libraries
library("limma")
library("data.table")

#### configs
# input
lmfit_object_path <- snakemake@input[["lmfit_object"]]
model_matrix_path <- snakemake@input[["model_matrix"]]
feature_annotation_path <- snakemake@input[["feature_annotation"]]

# output
contrast_result_path <- snakemake@output[["contrast_results"]]

# parameters
ova_var <- snakemake@wildcards[["ova_var"]]
feature_annotation_col <- base::make.names(snakemake@config[["feature_annotation"]][["column"]])[1]
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

#### identify groups for one-vs-all contrasts using column prefix & design matrix
groups <- colnames(design)[startsWith(colnames(design), ova_var)]

#### define contrasts
contrasts_all <- list()

for (gr in groups){
    # remove prefix
    gr_name <- gsub(ova_var, "", gr)
    # list of all groups without current group gr
    groups_wo_gr <- groups[groups != gr]
    # define contrast for one-vs-all analysis
    contrasts_all[[gr_name]] <- paste0(gr, " - (", paste(groups_wo_gr, collapse=" + "), ") / ", length(groups_wo_gr))
    print(contrasts_all[[gr_name]])
}

# generate contrast matrix based on contrast formulas and design matrix
contrast_matrix <- makeContrasts(contrasts=contrasts_all, levels = design)

#### fit contrasts
fit2 <- contrasts.fit(fit, contrast_matrix)

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
