
# visualize DEA results using volcano plots
rule volcanos:
    input:
        dea_results = os.path.join(result_path,'{analysis}','results.csv'),
    output:
        dea_volcanos = report(directory(os.path.join(result_path,'{analysis}','plots','volcano','{feature_list}')),
                              patterns=["{group}.png"],
                              caption="../report/volcano.rst",
                              category="{}_{}".format(config["project_name"], module_name),
                              subcategory="{analysis}",
                              labels={
                                  "name": "Volcano",
                                  "type": "{group}",
                                  "misc": "{feature_list}",
                              }),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/volcanos.yaml"
    log:
        os.path.join("logs","rules","volcanos_{analysis}_{feature_list}.log"),
    params:
        pCutoff = config["volcano"]["pCutoff"],
        FCcutoff = config["volcano"]["FCcutoff"],
        utils_path = workflow.source_path("../scripts/utils.R"),
    script:
        "../scripts/volcanos.R"
        
# visualize LFC of DEA results
rule lfc_heatmap:
    input:
        lambda wildcards: set() if wildcards.feature_list == 'FILTERED' else config["feature_lists"][wildcards.feature_list],
        dea_results = os.path.join(result_path,'{analysis}','results.csv'),
    output:
        dea_lfc_heatmap = #report( # report contains the following warning everywhere, hence excluded from report: <string>:191: (WARNING/2) Duplicate explicit target name: "filtered.png".
            os.path.join(result_path,'{analysis}','plots','heatmap','{feature_list}.png'),
                              # caption="../report/lfc_heatmap.rst",
                              # category="{}_{}".format(config["project_name"], module_name),
                              # subcategory="{analysis}",
                              # labels={
                              #     "name": "Heatmap",
                              #     "type": "Effect sizes",
                              #     "misc": "{feature_list}",
                              # }),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/heatmap.yaml"
    log:
        os.path.join("logs","rules","lfc_heatmap_{analysis}_{feature_list}.log"),
    params:
        adj_pval = config["filters"]["adj_pval"],
        lfc = config["filters"]["lfc"],
        ave_expr = config["filters"]["ave_expr"],
        utils_path = workflow.source_path("../scripts/utils.R"),
    script:
        "../scripts/heatmap.R"