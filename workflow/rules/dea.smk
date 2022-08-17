
# perform differential expression analysis
rule dea:
    input:
        get_paths
    output:
        dea_results = os.path.join(result_path,'{analysis}','DEA_results.csv'),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/limma.yaml"
    log:
        os.path.join("logs","rules","dea_{analysis}.log"),
    params:
        partition=config.get("partition"),
        feature_annotation = config["feature_annotation"],
        reference_levels = config["reference_levels"],
        formula = lambda w: annot_dict["{}".format(w.analysis)]["formula"],
        block_var = lambda w: annot_dict["{}".format(w.analysis)]["block_var"],
        comparisons = lambda w: annot_dict["{}".format(w.analysis)]["comparisons"],
        calcNormFactors_method = lambda w: annot_dict["{}".format(w.analysis)]["calcNormFactors_method"],
        voom = lambda w: annot_dict["{}".format(w.analysis)]["voom"],
        eBayes = lambda w: annot_dict["{}".format(w.analysis)]["eBayes"],
        limma_trend = lambda w: annot_dict["{}".format(w.analysis)]["limma_trend"],
    script:
        "../scripts/limma.R"

# aggregate results per analysis
rule aggregate:
    input:
        dea_results = os.path.join(result_path,'{analysis}','DEA_results.csv'),
    output:
        dea_stats = report(os.path.join(result_path,'{analysis}','DEA_stats.csv'), 
                                  caption="../report/dea_stats.rst", 
                                  category="{}_dea_limma".format(config["project_name"]), 
                                  subcategory="{analysis}"),
        dea_filtered_lfc = os.path.join(result_path,'{analysis}','DEA_FILTERED_LFC.csv'),
        dea_stats_plot = report(os.path.join(result_path,'{analysis}','plots','DEA_stats.png'), 
                                  caption="../report/dea_stats.rst", 
                                  category="{}_dea_limma".format(config["project_name"]), 
                                  subcategory="{analysis}"),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/ggplot.yaml"
    log:
        os.path.join("logs","rules","aggregate_{analysis}.log"),
    params:
        partition=config.get("partition"),
        adj_pval = config["filters"]["adj_pval"],
        lfc = config["filters"]["lfc"],
        ave_expr = config["filters"]["ave_expr"],
    script:
        "../scripts/aggregate.R"

