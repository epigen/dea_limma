#### load libraries & utility function 
library("limma")
library("edgeR")
library("statmod")
library("data.table")

# inputs
data_path <- snakemake@input[[1]]
metadata_path <- snakemake@input[[2]]
feature_annotation_path <- snakemake@input[["feature_annotation"]]

# outputs
dea_result_path <- snakemake@output[["dea_results"]]
lmfit_object_path <- snakemake@output[["lmfit_object"]]
model_matrix_path <- snakemake@output[["model_matrix"]]

# parameters
feature_annotation_col <- base::make.names(snakemake@params[["feature_annotation_col"]])[1] # make R compatible column name
reference_levels <- snakemake@params[["reference_levels"]] #list(treatment="UT", time="0h")
design <- formula(snakemake@params[["formula"]]) #formula("~ treatment")
block_var <- snakemake@params[["block_var"]] #0
comparisons <- strsplit(snakemake@params[["comparisons"]], "|", fixed=T)[[1]] #strsplit("treatment", "|", fixed=T)
calcNormFactors_method <- snakemake@params[["calcNormFactors_method"]] #"none"
voom_flag <- snakemake@params[["voom"]] #1
eBayes_flag <- snakemake@params[["eBayes"]] #1
limma_trend <- snakemake@params[["limma_trend"]] #0

result_dir <- dirname(dea_result_path)
# make directories if not exist
if (!dir.exists(result_dir)){
        dir.create(result_dir, recursive = TRUE)
    }

### load data
data <- data.frame(fread(file.path(data_path), header=TRUE), row.names=1)
stopifnot(isNumeric(data))
print("data")
print(dim(data))

### load metadata
metadata <- data.frame(fread(file.path(metadata_path), header=TRUE, na.strings=c("NA", "")), row.names=1)
rownames(metadata) <- make.names(rownames(metadata))
print("metadata")
print(dim(metadata))

### load feature annotation file (optional)
if (length(feature_annotation_path)!=0){
    feature_annot <- data.frame(fread(file.path(feature_annotation_path), header=TRUE), row.names=1)
    print("feature_annot")
    print(dim(feature_annot))
}

### prepare DEA

# subset metadata with used columns and non-NA rows
metadata_cols <- labels(terms(design)) #alternative: all.vars(design)
metadata_cols <- unique(unlist(lapply(metadata_cols, function(x) strsplit(x, ":", fixed = TRUE))))

if (block_var!=0){
    metadata_cols <- c(metadata_cols, block_var)
}
print("relevant metadata")
print(metadata_cols)

metadata <- metadata[,metadata_cols,drop=FALSE]
metadata <- na.omit(metadata)
print("metadata")
print(dim(metadata))

# keep intersection of data and metadata samples in same order
samples <- intersect(colnames(data), rownames(metadata))
data <- data[, samples]
metadata <- metadata[samples,,drop=FALSE]
stopifnot(all(colnames(data)==rownames(metadata)))
                                      
print("number of samples")
print(length(samples))

# relevel all categorical variables using the provided reference levels
for(col in names(reference_levels)){
    if(col %in% colnames(metadata)){
        metadata[[col]] <- as.factor(metadata[[col]])
        metadata[[col]] <- droplevels(metadata[[col]])
        metadata[[col]] <- relevel(metadata[[col]], ref=reference_levels[[col]])
    }
}

## MODEL
# create model matrix
model_matrix <- model.matrix(design, metadata)
# save model matrix
fwrite(as.data.frame(model_matrix), file=file.path(model_matrix_path), row.names=TRUE)

# check if the model represented by the design matrix has full rank ie linearly independent columns, which is required to make the model identifiable!
# new: more efficient and computationally stable compared to the svd() function, especially for large matrices, and it does not require rounding the singular values or checking for non-zero values.
# if(qr(model_matrix)$rank != ncol(model_matrix)){ # previous way of checking
if(!is.fullrank(model_matrix)){
    stop(paste0("Error: The design matrix representing your model does not have full rank, rendering the model not identifiable. The following columns are linearly dependent on previous columns: ", nonEstimable(model_matrix)))
}
# old: checks if all the singular values in the SVD decomposition of the model_matrix are non-zero, indicating that the matrix is full rank and invertible.
# stopifnot(all(round(svd(model_matrix)$d, 6) != 0))

       
## PREPARE OBJECTS
# calculate Normalization Factors (optional)
if (calcNormFactors_method!=0){
    # create dge object
    dge <- DGEList(data, samples=metadata, genes=rownames(data))
    dge <- calcNormFactors(dge, method=calcNormFactors_method)
    
    # voom (optional)
    if (voom_flag!=0){
        pdf(file=file.path(result_dir,"mean_var_trend_voom.pdf"))
        v <- voom(dge, model_matrix, plot=TRUE)
        x <- dev.off()
    }else{
        # directly taken from "Differential expression: limma-trend" in the limma userguide
        # http://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf
        v <- edgeR::cpm(dge, log = TRUE, prior.count = 3)
    }
}else{ # ie if calcNormFactors_method==0 -> data is already log-normalized
     v <- data
}


### perform DEA

# fit linear model with blocking (optional) or without
if(block_var!=0){
    block <- metadata[[block_var]]
    
    print("duplicateCorrelation")
    corr_fit <- duplicateCorrelation(v, model_matrix, block=block)
    print(corr_fit$consensus)
    cons_correlation <- corr_fit$consensus
    
    # voom (optional)
    if (voom_flag!=0){
        # now it is potentially beneficial to voom again (slower)
        print("vooming")
        pdf(file=file.path(result_dir,"mean_var_trend_voom2.pdf"))
        v <- voom(dge, model_matrix, block=block, correlation=cons_correlation, plot=TRUE)
        x <- dev.off()

        # now it is potentially beneficial to calculate correlation again
        print("duplicateCorrelation")
        corr_fit <- duplicateCorrelation(v, model_matrix, block=block)
        print(corr_fit$consensus)
        cons_correlation <- corr_fit$consensus
    }
    
    print("fitting")
    lmfit = lmFit(v, model_matrix, block=block, correlation=cons_correlation)
} else {
    lmfit <- lmFit(v, model_matrix)
}
                                      
# save fitted model object for manual downstream analyses e.g., contrasts
saveRDS(lmfit, file = file.path(lmfit_object_path))
# save the used expression matrix input object
saveRDS(v, file = file.path(result_dir, "expression_matrix.rds"))

# plot mean-variance trend of fitted model
pdf(file=file.path(result_dir,"mean_var_trend_fitted.pdf"))
plotSA(lmfit)
x <- dev.off()

# eBayes (optional) without or with limma-trend (optional)
if (eBayes_flag!=0){
    if(limma_trend==0){
        lmfit <- eBayes(lmfit, robust=TRUE, trend=FALSE)
    }else{
        lmfit <- eBayes(lmfit, robust=TRUE, trend=TRUE)
    }
}else{ # Warning: not recommended by limma author Gordon Smyth (https://support.bioconductor.org/p/35174/)
    # determine ordinary t statistic, p-value and B-statistic (lods) according to
    # limma userguide chapter 13.2 Fitted Model Objects: http://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf
    # bioconductor forum: https://support.bioconductor.org/p/113833/
    # bioconductor forum: https://support.bioconductor.org/p/35159/
    # note: B-statistic determined using eBayes, assuming it will not be used downstream
    
    lmfit$t <- lmfit$coef/lmfit$stdev.unscaled/lmfit$sigma
    lmfit$p.value <- 2 * pt(-abs(lmfit$t), df = lmfit$df.residual)
    lmfit$lods <- eBayes(lmfit, robust=TRUE, trend=FALSE)$lods
}

# topTable to extract coefficients
dea_results <- data.frame()

for(coefx in colnames(coef(lmfit))){
    if(any(unlist(lapply(comparisons, function(x) grepl(x, coefx))))){
        print(coefx)
        
        tmp_res <- topTable(lmfit, coef=coefx, number=nrow(data), sort.by="P")
        tmp_res$feature <- rownames(tmp_res)
        tmp_res <- tmp_res[,c(ncol(tmp_res),1:(ncol(tmp_res)-1))]
        rownames(tmp_res) <- NULL
        tmp_res$group <- gsub(":", "__", coefx) # replace colon with double-underscore
        
        if (length(feature_annotation_path)!=0){
            tmp_res$feature_name <- feature_annot[tmp_res$feature, feature_annotation_col]
        }
    
        if(dim(dea_results)[1]==0){
            dea_results <- tmp_res
        }else{
            dea_results <- rbind(dea_results, tmp_res)
        }
    }
}
                         
# remove rows with adj.P.Val=NA
dea_results <- dea_results[!is.na(dea_results$adj.P.Val),]

### save DEA results
fwrite(as.data.frame(dea_results), file=file.path(dea_result_path), row.names=FALSE)
