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
ova_var <- snakemake@params[["original_ova_var"]]
feature_annotation_col <- base::make.names(snakemake@params[["feature_annotation_col"]])[1]
eBayes_flag <- snakemake@params[["eBayes"]]
limma_trend <- snakemake@params[["limma_trend"]]
formula <- snakemake@params[["formula"]]

#### assert that ova_var is not an interaction term
if (grepl(":", ova_var, fixed=TRUE)) {
    stop(paste0("One-vs-all contrasts are not supported for interaction terms. '", ova_var, "' is an interaction term."))
}

#### load data (model, design, feature_annotation and metadata)
fit <- readRDS(lmfit_object_path)
design <- data.frame(fread(file.path(model_matrix_path), header=TRUE), row.names=1)
meta <- data.frame(fread(file.path(metadata_path), header=TRUE), row.names=1)

# fit still has : in the column names since it was read from a RDS, but the design matrix replaced them
# since it is a dataframe (read.csv would do so as well) --> adapt the fit object to match
colnames(fit$coefficients) <- make.names(colnames(fit$coefficients))

### load feature annotation file (optional)
if (length(feature_annotation_path)!=0){
    feature_annot <- data.frame(fread(file.path(feature_annotation_path), header=TRUE), row.names=1)
    print("feature_annot")
    print(dim(feature_annot))
}

### determine flags

# find out if there is a term that was modelled as a means model and if yes, which metadata column that was
# i.e., which is the first non-interaction term in the formula
is_means_model <- attributes(terms(formula(formula)))$intercept == 0
all_terms <- attributes(terms(formula(formula)))$term.labels
means_model_term <- FALSE
if (is_means_model){    
    for (i in 1:length(all_terms)){
        term <- all_terms[i]
        if (grepl(":", term, fixed=TRUE)){
            # interaction term, skip
            next
        }
        means_model_term <- gsub(" ", "", term)
        break
    }
}
print(paste0("Means model term: ", means_model_term))
# is the current OvA terms is the one that's modelled as a means model
ova_is_means_model <- ova_var == means_model_term

#### identify groups for one-vs-all contrasts using the metadata table, column prefix & design matrix
# check if ova_var is numeric, and if yes, cause error since OvA does not make sense for continuous variables
# ok to kill the process, since the user should not get to make this wrong choice
if (is.numeric(meta[, ova_var])){
    stop(paste0("One-vs-all contrasts do not make sense for numeric (non-factor) variable '", ova_var, "'"))
}

group_names <- unique(meta[,ova_var])
group_cols <- base::make.names(paste0(ova_var, group_names))

#### define contrasts
contrasts_all <- list()
for (i in seq(length(group_names))){
    group_name <- group_names[[i]]
    group_col <- group_cols[[i]]
    group_cols_wo_gr <- group_cols[group_cols != group_col]

    # see https://github.com/epigen/dea_limma/issues/34 for a more detailed explanation of what is happening here
    if (ova_is_means_model){
        # for a means model, simply subtract the mean of all other groups from the current group
        contrasts_all[[group_name]] <- paste0(
            group_col, " - (", paste(group_cols_wo_gr, collapse=" + "), ") / ", length(group_cols_wo_gr)
        )
    
    } else {
        # each actual effect of the levels contains the intercept, so in the when subtracting the average of all other
        # groups from the current group the intercept cancels out. Thus, we compare relative effects (relative to the 
        # reference group), where the relative effect for the reference group itself is 0 (in the contrast formula just
        # falls away, but is counted for the denominator in the mean).
    
        reference_group <- setdiff(group_cols, colnames(design))
        if (length(reference_group) != 1) {
            # something went wrong
            stop("There should be exactly one reference group in the design matrix.")
        } else {
            reference_group <- reference_group[1]
            print(paste0("Reference group: ", reference_group))
        }

        self_relative_effect = ifelse(group_col == reference_group, "", group_col)
        
        # collect the effects of all the other groups and sum then up, to then take the average
        other_group_effects <- c()
        for (other_group_col in group_cols_wo_gr){
            if (other_group_col == reference_group) {
                # if the other group is the reference group, it does not contribute to the contrast formula
                # so we just skip it (don't add empty string to avoid `+  +`)
                next
            }
            other_group_effects <- c(other_group_effects, other_group_col)
        }
        
        if (length(other_group_effects) == 0){
            # if there are only two levels and one is the reference, only one term remains for the contrast formula
            contrasts_all[[group_name]] <- self_relative_effect
        } else {
            other_group_average <- paste0(
                "(", paste(other_group_effects, collapse=" + "), ") / ", length(group_cols_wo_gr)
            )
            contrasts_all[[group_name]] <- paste0(self_relative_effect, " - ", other_group_average)        
        }
    }
    
    print(paste0("Contrast for group '", group_name, "': "))
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
    if (length(feature_annotation_path)!=0){
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
contrast_result$group <- gsub("\\.", "__", contrast_result$group)
fwrite(as.data.frame(contrast_result), file=file.path(contrast_result_path), row.names=FALSE)