#### load libraries & utility function 
library(EnhancedVolcano, quietly=TRUE)
library(patchwork, quietly=TRUE)
library(ggplot2)

# source utility functions
source("workflow/scripts/utils.R")

# inputs
dea_result_path <- snakemake@input[["dea_results"]] #"/nobackup/lab_bock/projects/macroIC/results/CC001/dea_limma/CC001_thp1_filtered/DEA_results.csv" 

# outputs
volcano_plot_path <- snakemake@output[["dea_volcanos"]] #"/nobackup/lab_bock/projects/macroIC/results/CC001/dea_limma/CC001_thp1_filtered/plots/DEA_volcanos.png"

# parameters
pCutoff <- snakemake@params[["pCutoff"]] # 0.05
FCcutoff <- snakemake@params[["FCcutoff"]] # 2


# plot specifications
width <- 4
height <- 5

### load DEA results
dea_results <- read.csv(file=file.path(dea_result_path))

n_col <- min(10, length(unique(dea_results$group)))
width_panel <- n_col * width + 1
height_panel <- height * ceiling(length(unique(dea_results$group))/n_col)

### Visualize DEA results using Volcano plots

volcano_plots <- list()

for (group in unique(dea_results$group)){
    toptable <- dea_results[dea_results$group==group,]
    
    if("feature_name" %in% colnames(toptable)){
        lab <- toptable$feature_name
    }else{
        lab <- toptable$feature
    }
    
    x <- "logFC"
    y <- "adj.P.Val"

    volcano_plots[[group]] <- EnhancedVolcano(toptable = toptable,
                    lab = lab,
                    x = x,
                    y = y,
                    selectLab = NULL,
                    xlim = c(min(toptable[[x]], na.rm = TRUE) - 1, max(toptable[[x]], na.rm = TRUE) + 1),
                    ylim = c(0, max(-log10(toptable[[y]]), na.rm = TRUE) + 5),
                    xlab = bquote(~log[2] ~ "fold change"),
                    ylab = bquote(~-log[10] ~ "adjusted p-value"),
                    axisLabSize = 12,
                    title = paste0(group),
                    subtitle = '', # default: bquote(italic(EnhancedVolcano))
                    caption = paste0("variables:",nrow(toptable),"; log2FC>",FCcutoff,"; adj.p-val<",pCutoff),
                    titleLabSize = 14,
                    subtitleLabSize = 0,
                    captionLabSize = 6,
                    pCutoff = pCutoff, #default: 0.05
                    pCutoffCol = y,
                    FCcutoff = FCcutoff, # default:1
                    cutoffLineType = "longdash",
                    cutoffLineCol = "black",
                    cutoffLineWidth = 0.4,
                    pointSize = 1, # default: 2
                    labSize = 3, #default: 5
                    labCol = "black",
                    labFace = "plain",
                    boxedLabels = TRUE, #default: FALSE
                    parseLabels = FALSE,
                    shape = 19,
                    shapeCustom = NULL,
                    col = c("grey30", "forestgreen", "royalblue", "red2"),
                    colCustom = NULL,
                    colAlpha = 1/2,
                    colGradient = NULL,
                    colGradientBreaks = c(pCutoff, 1),
                    colGradientLabels = c("0", "1.0"),
                    colGradientLimits = c(0, 1),
                    legendLabels = c("NS", expression(~ log[2] ~ FC), "adj. p-value", 'both'),
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
                    drawConnectors = TRUE, #default: FALSE
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

}

volcano_plots_panel <- wrap_plots(volcano_plots, ncol = n_col, guides = "collect")

# save plot (plotting in jupyterLab does not work, but ggsave works) https://github.com/slowkow/ggrepel/issues/113
# options(repr.plot.width=width_panel, repr.plot.height=height_panel)
# print(volcano_plots_panel)

set.seed(42)
ggsave_new(filename = "DEA_volcanos", 
           results_path=dirname(volcano_plot_path), 
           plot=volcano_plots_panel, 
           width=width_panel, 
           height=height_panel)
