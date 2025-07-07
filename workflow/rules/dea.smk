
# perform differential expression analysis
rule dea:
    input:
        get_paths,
        feature_annotation = config["feature_annotation"]["path"] if config["feature_annotation"]["path"]!="" else [],
    output:
        dea_results = os.path.join(result_path,'{analysis}','results.csv'),
        lmfit_object = os.path.join(result_path,'{analysis}','lmfit_object.rds'),
        model_matrix = os.path.join(result_path,'{analysis}','model_matrix.csv'),
    wildcard_constraints:
        analysis = "(?!.*_OvA_).*"
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/limma.yaml"
    log:
        os.path.join("logs","rules","dea_{analysis}.log"),
    params:
        feature_annotation_col = config["feature_annotation"]["column"],
        reference_levels = config["reference_levels"],
        formula = lambda w: annot_dict["{}".format(w.analysis)]["formula"],
        block_var = lambda w: annot_dict["{}".format(w.analysis)]["block_var"],
        # protect against empty comparison variable
        comparisons = lambda w: annot_dict[str(w.analysis)]["comparisons"] if isinstance(annot_dict[str(w.analysis)]["comparisons"], str) else "",
        calcNormFactors_method = lambda w: annot_dict["{}".format(w.analysis)]["calcNormFactors_method"],
        voom = lambda w: annot_dict["{}".format(w.analysis)]["voom"],
        eBayes = lambda w: annot_dict["{}".format(w.analysis)]["eBayes"],
        limma_trend = lambda w: annot_dict["{}".format(w.analysis)]["limma_trend"],
    script:
        "../scripts/limma.R"

# performs a one-vs-all (OvA) differential analysis using contrasts based on the previous model
rule one_vs_all_contrasts:
    input:
        metadata = get_metadata_path,
        lmfit_object = os.path.join(result_path,'{analysis}','lmfit_object.rds'),
        model_matrix = os.path.join(result_path,'{analysis}','model_matrix.csv'),
        feature_annotation = config["feature_annotation"]["path"] if config["feature_annotation"]["path"]!="" else [],
    output:
        contrast_results = os.path.join(result_path,'{analysis}_OvA_{ova_var}','results.csv'),
        contrast_object = os.path.join(result_path,'{analysis}_OvA_{ova_var}','lmfit_object.rds'),
        contrast_matrix = os.path.join(result_path,'{analysis}_OvA_{ova_var}','model_matrix.csv'),
    params:
        eBayes = lambda w: annot_dict["{}".format(w.analysis)]["eBayes"],
        limma_trend = lambda w: annot_dict["{}".format(w.analysis)]["limma_trend"],
        feature_annotation_col = config["feature_annotation"]["column"],
        formula = lambda w: annot_dict["{}".format(w.analysis)]["formula"],
        reference_levels = config["reference_levels"],
        original_ova_var = lambda w: ova_analyses[f"{w.analysis}_OvA_{w.ova_var}"]
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/limma.yaml"
    log:
        os.path.join("logs","rules","one_vs_all_contrasts_{analysis}_{ova_var}.log"),
    script:
        "../scripts/one_vs_all_contrasts.R"

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
        utils_path=workflow.source_path("../scripts/utils.R"),
    script:
        "../scripts/aggregate.R"

rule ova_stats_plot:
    input:
        dea_results = get_input_ova_stats_plot,
    output:
        dea_stats_OvA = report(os.path.join(result_path,'{analysis}','stats_OvA.csv'), 
                                  caption="../report/dea_stats_OvA.rst", 
                                  category="{}_{}".format(config["project_name"], module_name),
                                  subcategory="{analysis}",
                                  labels={
                                      "name": "Statistics",
                                      "type": "table",
                                      "misc": "CSV",
                                  }),
        dea_stats_plot_OvA = report(os.path.join(result_path,'{analysis}','plots','stats_OvA.png'),
                                  caption="../report/dea_stats_OvA.rst", 
                                  category="{}_{}".format(config["project_name"], module_name),
                                  subcategory="{analysis}",
                                  labels={
                                      "name": "Statistics",
                                      "type": "Bar plot",
                                      "misc": "PNG",
                                  }),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/ggplot.yaml"
    log:
        os.path.join("logs","rules","ova_stats_plot_{analysis}.log"),
    params:
        score_formula = config["score_formula"],
        adj_pval = config["filters"]["adj_pval"],
        lfc = config["filters"]["lfc"],
        ave_expr = config["filters"]["ave_expr"],
        utils_path=workflow.source_path("../scripts/utils.R"),
    script:
        "../scripts/ova_stats_plot.R"

# shadow rule to enable downstream processing
# requires to know that the file will exist in that exact location, otherwise MissingInputException error
rule fetch_file:
    input:
        dea_stats = os.path.join(result_path,'{analysis}','stats.csv'),
    output:
        feature_list = update(os.path.join(result_path,'{analysis}','feature_lists',"{group}_{type,(up_features.txt|up_features_annot.txt|down_features.txt|down_features_annot.txt|featureScores.csv|featureScores_annot.csv)}")),
    resources:
        mem_mb="1000",
    shell:
        """
        # only if the file already exists
        if [ -f {output.feature_list} ]; then \
            touch {output.feature_list}; \
        fi
        """
