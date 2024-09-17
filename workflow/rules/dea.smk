
# perform differential expression analysis
rule dea:
    input:
        get_paths
    output:
        dea_results = os.path.join(result_path,'{analysis}','results.csv'),
        lmfit_object = os.path.join(result_path,'{analysis}','lmfit_object.rds'),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/limma.yaml"
    log:
        os.path.join("logs","rules","dea_{analysis}.log"),
    params:
        feature_annotation = config["feature_annotation"],
        reference_levels = config["reference_levels"],
        formula = lambda w: annot_dict["{}".format(w.analysis)]["formula"],
        formula_batch = lambda w: annot_dict["{}".format(w.analysis)]["formula_batch"],
        block_var = lambda w: annot_dict["{}".format(w.analysis)]["block_var"],
        batch_var = lambda w: annot_dict["{}".format(w.analysis)]["batch_var"],
        comparisons = lambda w: annot_dict["{}".format(w.analysis)]["comparisons"],
        calcNormFactors_method = lambda w: annot_dict["{}".format(w.analysis)]["calcNormFactors_method"],
        filterByExpr = lambda w: annot_dict["{}".format(w.analysis)]["filterByExpr"],
        voom = lambda w: annot_dict["{}".format(w.analysis)]["voom"],
        eBayes = lambda w: annot_dict["{}".format(w.analysis)]["eBayes"],
        limma_trend = lambda w: annot_dict["{}".format(w.analysis)]["limma_trend"],
    script:
        "../scripts/limma.R"

# aggregate results per analysis
checkpoint aggregate:
    input:
        dea_results = os.path.join(result_path,'{analysis}','results.csv'),
    output:
        all_features = os.path.join(result_path,'{analysis}','feature_lists','ALL_features.txt'),
        feature_lists = directory(os.path.join(result_path,'{analysis}','feature_lists')),
        dea_stats = report(os.path.join(result_path,'{analysis}','stats.csv'), 
                                  caption="../report/dea_stats.rst", 
                                  category="{}_{}".format(config["project_name"], module_name),
                                  subcategory="{analysis}",
                                  labels={
                                      "name": "Statistics",
                                      "type": "table",
                                      "misc": "CSV",
                                  }),
        dea_stats_plot = report(os.path.join(result_path,'{analysis}','plots','stats.png'), 
                                  caption="../report/dea_stats.rst", 
                                  category="{}_{}".format(config["project_name"], module_name),
                                  subcategory="{analysis}",
                                  labels={
                                      "name": "Statistics",
                                      "type": "Bar plot",
                                      "misc": "PNG",
                                  }),
        dea_pvalue_plot = report(directory(os.path.join(result_path,'{analysis}','plots','pvalue_distribution')), 
                                  patterns=["{group}.png"],
                                 caption="../report/pvalue_distribution.rst", 
                                  category="{}_{}".format(config["project_name"], module_name),
                                  subcategory="{analysis}",
                                  labels={
                                      "name": "P-value distribution",
                                      "type": "{group}",
                                      "misc": "PNG",
                                  }),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/ggplot.yaml"
    log:
        os.path.join("logs","rules","aggregate_{analysis}.log"),
    params:
        score_formula = config["score_formula"],
        adj_pval = config["filters"]["adj_pval"],
        lfc = config["filters"]["lfc"],
        ave_expr = config["filters"]["ave_expr"],
    script:
        "../scripts/aggregate.R"

# intermediate rule required for checkpoint usage with wildcards to aggregate all result feature lists
rule gather_feature_lists:
    input:
        get_feature_lists,
    output:
        os.path.join(result_path, '{analysis}', 'result_feature_lists.txt'),
    resources:
        mem_mb="1000",
    threads: 1
    log:
        os.path.join("logs","rules","gather_feature_lists_{analysis}.log"),
    shell:
        """
        echo '{input}' | tr ' ' '\n' > {output}
        """