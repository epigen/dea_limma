
# visualize DEA results using volcano plots
rule volcanos:
    input:
        dea_results = os.path.join(result_path,'{analysis}','DEA_results.csv'),
    output:
        dea_volcanos = report(os.path.join(result_path,'{analysis}','plots','DEA_volcanos_{feature_list}_{pval_type}.png'),
                              caption="../report/volcano.rst",
                              category="{}_{}".format(config["project_name"], module_name),
                              subcategory="{analysis}",
                              labels={
                                  "name": "Volcano plot",
                                  "type": "{pval_type}",
                                  "misc": "{feature_list}",
                              }),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/volcanos.yaml"
    log:
        os.path.join("logs","rules","volcanos_{analysis}_{feature_list}_{pval_type}.log"),
    params:
        partition=config.get("partition"),
        pCutoff = config["volcano"]["pCutoff"],
        FCcutoff = config["volcano"]["FCcutoff"],
    script:
        "../scripts/volcanos.R"
        
# visualize LFC of DEA results
rule lfc_heatmap:
    input:
        dea_results = os.path.join(result_path,'{analysis}','DEA_results.csv'),
#         dea_lfc = os.path.join(result_path,'{analysis}','DEA_LFC.csv'),
    output:
        dea_lfc_heatmap = report(os.path.join(result_path,'{analysis}','plots','DEA_LFC_heatmap_{feature_list}.png'),
                              caption="../report/lfc_heatmap.rst",
                              category="{}_{}".format(config["project_name"], module_name),
                              subcategory="{analysis}",
                              labels={
                                  "name": "Heatmap",
                                  "type": "effect sizes",
                                  "misc": "{feature_list}",
                              }),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/heatmap.yaml"
    log:
        os.path.join("logs","rules","lfc_heatmap_{analysis}_{feature_list}.log"),
    params:
        partition=config.get("partition"),
        adj_pval = config["filters"]["adj_pval"],
        lfc = config["filters"]["lfc"],
        ave_expr = config["filters"]["ave_expr"],
    script:
        "../scripts/heatmap.R"