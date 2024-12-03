
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
        dea_results = os.path.join(result_path,'{analysis}','results.csv'),
    output:
        all_features = os.path.join(result_path,'{analysis}','feature_lists','ALL_features.txt'),
        all_features_annot = os.path.join(result_path,'{analysis}','feature_lists','ALL_features_annot.txt') if config["feature_annotation"]["path"]!="" else [],
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
        group = "ALL",
        score_formula = config["score_formula"],
        adj_pval = config["filters"]["adj_pval"],
        lfc = config["filters"]["lfc"],
        ave_expr = config["filters"]["ave_expr"],
    script:
        "../scripts/aggregate.R"

# generate feature lists per group
# not part of target rule all and only used when the outputs are required by a subsequent module e.g., enrichment_analysis
rule feature_lists:
    input:
        dea_results = os.path.join(result_path,'{analysis}','results.csv'),
        dea_stats = os.path.join(result_path,'{analysis}','stats.csv'), # this ensures rule order (after rule aggregate) to avoid file writing clashes
    output:
        features_up = os.path.join(result_path,'{analysis}','feature_lists',"{group}_up_features.txt"),
        features_up_annot = os.path.join(result_path,'{analysis}','feature_lists',"{group}_up_features_annot.txt") if config["feature_annotation"]["path"]!="" else [],
        features_down = os.path.join(result_path,'{analysis}','feature_lists',"{group}_down_features.txt"),
        features_down_annot = os.path.join(result_path,'{analysis}','feature_lists',"{group}_down_features_annot.txt") if config["feature_annotation"]["path"]!="" else [],
        features_scores = os.path.join(result_path,'{analysis}','feature_lists',"{group}_featureScores.csv") if config["score_formula"]!="" else [],
        features_scores_annot = os.path.join(result_path,'{analysis}','feature_lists',"{group}_featureScores_annot.csv") if (config["score_formula"]!="" and config["feature_annotation"]["path"]!="") else [],
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/ggplot.yaml"
    log:
        os.path.join("logs","rules","feature_lists_{analysis}_{group}.log"),
    params:
        group = lambda w: "{}".format(w.group),
        score_formula = config["score_formula"],
        adj_pval = config["filters"]["adj_pval"],
        lfc = config["filters"]["lfc"],
        ave_expr = config["filters"]["ave_expr"],
    script:
        "../scripts/aggregate.R"