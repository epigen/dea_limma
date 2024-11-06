#### load libraries & utility function 
library("EnhancedVolcano", quietly=TRUE)
# library("patchwork", quietly=TRUE)
library("ggplot2")
library("data.table")

# source utility functions
# source("workflow/scripts/utils.R")
snakemake@source("utils.R")

# inputs
dea_result_path <- snakemake@input[["dea_results"]]

# outputs
volcano_plot_path <- snakemake@output[["dea_volcanos"]]

# parameters
pCutoff <- snakemake@params[["pCutoff"]] # 0.05
FCcutoff <- snakemake@params[["FCcutoff"]] # 2
pval_type <- snakemake@wildcards[["pval_type"]] # 'adjp'
feature_list_name <- snakemake@wildcards[["feature_list"]]

# plot specifications
width <- 5
height <- 4

### load DEA results
# dea_results <- read.csv(file=file.path(dea_result_path))
dea_results <- data.frame(fread(file.path(dea_result_path), header=TRUE))

# load feature list if not ALL
if (feature_list_name!="ALL"){
    feature_list_path <- snakemake@config[["feature_lists"]][[feature_list_name]]
    feature_list <- scan(file.path(feature_list_path), character())
}

# replace (adjusted) p-value of 0 with min(adj.P.Val>0) * 10ˆ-1 for plotting
dea_results[dea_results$adj.P.Val==0,"adj.P.Val"] <- min(dea_results$adj.P.Val[dea_results$adj.P.Val>0])*10^-1
dea_results[dea_results$P.Value==0,"P.Value"] <- min(dea_results$P.Value[dea_results$P.Value>0])*10^-1

for (pval_type in c("adj.P.Val", "P.Value")){
    for (group in unique(dea_results$group)){
        
        toptable <- dea_results[dea_results$group==group,]

        # set volcano parameters
        x <- "logFC"
        selectLab <- NULL
        colCustom <- NULL
        colAlpha <- 1/2

        # if feature list provided then color features of interest
        if (feature_list_name!="ALL"){
            if("feature_name" %in% colnames(toptable)){
                toptable$feature_list <- ifelse(toptable$feature_name %in% feature_list, TRUE, FALSE)
                # sort results by feature list so they are on top in the plot
                toptable <- toptable[order(toptable$feature_list),]
                # highlight features of interest
                keyvals.alpha <- ifelse(toptable$feature_name %in% feature_list, 0.5, 0.25)
                keyvals.col <- ifelse(toptable$feature_name %in% feature_list, "red", "grey")
            }else{
                toptable$feature_list <- ifelse(toptable$feature %in% feature_list, TRUE, FALSE)
                # sort results by feature list so they are on top in the plot
                toptable <- toptable[order(toptable$feature_list),]
                # highlight features of interest
                keyvals.alpha <- ifelse(toptable$feature %in% feature_list, 0.5, 0.25)
                keyvals.col <- ifelse(toptable$feature %in% feature_list, "red", "grey")
            }

            # name features of interest
            names(keyvals.col)[keyvals.col == "red"] <- feature_list_name
            names(keyvals.col)[keyvals.col == "grey"] <- "other features"

            # set volcano parameters
            selectLab <- feature_list
            colCustom <- keyvals.col
            colAlpha <- keyvals.alpha
        }

		type_label <- if(pval_type=='adj.P.Val') "adj. p-value" else "raw p-value"
        volcano_plot <- EnhancedVolcano(toptable = toptable,
                        lab = if (feature_list_name=="ALL") NA else if("feature_name" %in% colnames(toptable)) toptable$feature_name else toptable$feature,
                        x = x,
                        y = pval_type,
                        selectLab = if (feature_list_name=="ALL") NULL else selectLab,
                        xlim = c(min(toptable[[x]], na.rm = TRUE) - 1, max(toptable[[x]], na.rm = TRUE) + 1),
                        ylim = c(0, max(-log10(toptable[[pval_type]]), na.rm = TRUE) + 5),
                        xlab = bquote(~log[2] ~ "fold change"),
                        ylab = if(pval_type=='adj.P.Val') bquote(~-log[10] ~ "adjusted p-value") else bquote(~-log[10] ~ "raw p-value"),
                        axisLabSize = 12,
                        title = addline_format(group),
                        subtitle = '', # default: bquote(italic(EnhancedVolcano))
                        caption = paste0("variables:",nrow(toptable),"; |avg.log2FC|>",FCcutoff,"; adj.p-val<",pCutoff),
                        titleLabSize = 14,
                        subtitleLabSize = 0,
                        captionLabSize = 6,
                        pCutoff = pCutoff, #default: 0.05
                        pCutoffCol = pval_type, #"adj.P.Val",
                        FCcutoff = FCcutoff, # default:1
                        cutoffLineType = "longdash",
                        cutoffLineCol = "black",
                        cutoffLineWidth = 0.4,
                        pointSize = 0.5, # default: 2
                        labSize = 2, #default: 5
                        labCol = "black",
                        labFace = "plain",
                        boxedLabels = TRUE, #default: FALSE
                        parseLabels = FALSE,
                        shape = 19,
                        shapeCustom = NULL,
                        col = c("grey30", "forestgreen", "royalblue", "red2"),
                        colCustom = colCustom,
                        colAlpha = colAlpha,
                        colGradient = NULL,
                        colGradientBreaks = c(pCutoff, 1),
                        colGradientLabels = c("0", "1.0"),
                        colGradientLimits = c(0, 1),
                        legendLabels = c("NS", expression(~ log[2] ~ FC), type_label, 'both'),
                        legendPosition = "right", #default: "top"
                        legendLabSize = 14,
                        legendIconSize = 4,
                        legendDropLevels = TRUE,
                        encircle = NULL,
                        encircleCol = "black",
                        encircleFill = "pink",
                        encircleAlpha = 3/4,
                        encircleSize = 2.5,
                        shade = NULL,
                        shadeFill = "grey",
                        shadeAlpha = 1/2,
                        shadeSize = 0.01,
                        shadeBins = 2,
                        drawConnectors = if (feature_list_name=="ALL") FALSE else TRUE, #default: FALSE
                        widthConnectors = 0.1, # default: 0.5
                        typeConnectors = "closed",
                        endsConnectors = "first",
                        lengthConnectors = unit(0.01, "npc"),
                        colConnectors = "grey10",
                        max.overlaps = 0,
                        maxoverlapsConnectors = 10, # default: NULL
                        min.segment.length = 0,
                        directionConnectors = "both",
                        arrowheads = FALSE, # default: TRUE
                        hline = NULL,
                        hlineType = "longdash",
                        hlineCol = "black",
                        hlineWidth = 0.4,
                        vline = NULL,
                        vlineType = "longdash",
                        vlineCol = "black",
                        vlineWidth = 0.4,
                        gridlines.major = TRUE,
                        gridlines.minor = TRUE,
                        border = "partial",
                        borderWidth = 0.8,
                        borderColour = "black",
                        raster = FALSE
                       ) + custom_theme + theme(legend.title = element_blank())
        
        set.seed(42)
        ggsave_new(filename = if(pval_type=='adj.P.Val') paste0(group,"_adjp") else paste0(group,"_rawp"),
                   results_path=volcano_plot_path, 
                   plot=volcano_plot, 
                   width=width, 
                   height=height)    
    }
}
