#### load libraries & utility function 
library(limma)
library(edgeR)
library(statmod)

# inputs
data_path <- snakemake@input[[1]] #"/nobackup/lab_bock/projects/macroIC/results/CC001/counts/thp1_filtered.csv"
metadata_path <- snakemake@input[[2]] #"/nobackup/lab_bock/projects/macroIC/results/CC001/counts/thp1_annotation.csv"

# outputs
dea_result_path <- snakemake@output[["dea_results"]] #"/nobackup/lab_bock/projects/macroIC/results/CC001/dea_limma/CC001_thp1_filtered/DEA_results.csv" 

# parameters
feature_annotation <- snakemake@params[["feature_annotation"]] #list(path="/nobackup/lab_bock/projects/macroIC/results/CC001/counts/gene_annotation.csv", column="external_gene_name")
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
data <- read.csv(file=data_path, row.names=1)
stopifnot(isNumeric(data))
print("data")
print(dim(data))

### load metadata
metadata <- read.csv(file=metadata_path, row.names=1, na.strings=c("NA",""))
rownames(metadata) <- make.names(rownames(metadata))
print("metadata")
print(dim(metadata))

### load feature annotation file (optional)
if (feature_annotation[["path"]]!=""){
    feature_annot <- read.csv(file=feature_annotation[["path"]], row.names=1)
    print("feature_annot")
    print(dim(feature_annot))
}

### prepare DEA

# subset metadata with used columns and non-NA rows
metadata_cols <- labels(terms(design)) #alternative all.vars(design)
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

# relevel all categorical variables using the provided reference levels
for(col in names(reference_levels)){
    if(col %in% colnames(metadata)){
        metadata[[col]] <- as.factor(metadata[[col]])
        metadata[[col]] <- droplevels(metadata[[col]])
        metadata[[col]] <- relevel(metadata[[col]], ref=reference_levels[[col]])
    }
}


# calculate Normalization Factors (optional)
if (calcNormFactors_method!=0){
    # create dge object
    dge <- DGEList(data, samples=metadata, genes=rownames(data))
    dge <- calcNormFactors(dge, method=calcNormFactors_method)
}else{
    dge <- data
}

# create model matrix
model_matrix <- model.matrix(design, metadata)
# save model matrix
write.csv(model_matrix, file=file.path(result_dir,"model_matrix.csv"))

# check if the model is not over-determined using SVD TODO: provide source or explanation
stopifnot(all(round(svd(model_matrix)$d, 6) != 0))

# voom (optional)
if (voom_flag!=0){
    pdf(file=file.path(result_dir,"mean_var_trend_voom.pdf"))
    v <- voom(dge, model_matrix, plot=TRUE)
    x <- dev.off()
}else{
    v <- dge
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
        
        if (feature_annotation[["path"]]!=""){
            tmp_res$feature_name <- feature_annot[tmp_res$feature, feature_annotation[["column"]]]
        }
    
        if(dim(dea_results)[1]==0){
            dea_results <- tmp_res
        }else{
            dea_results <- rbind(dea_results, tmp_res)
        }
    }
}

# rename "genes" column to "feature"
# colnames(dea_results)[colnames(dea_results)=="genes"] <- "feature"

### save DEA results
write.csv(dea_results, file=file.path(dea_result_path), row.names=FALSE)

