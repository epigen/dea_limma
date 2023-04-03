
# visualize DEA results using volcano plots
rule volcanos:
    input:
        dea_results = os.path.join(result_path,'{analysis}','DEA_results.csv'),
    output:
        dea_volcanos = report(os.path.join(result_path,'{analysis}','plots','DEA_volcanos_{feature_list}_{pval_type}.png'),
                              caption="../report/volcano.rst",
                              category="{}_dea_limma".format(config["project_name"]),
                              subcategory="{analysis}"),
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
        dea_filtered_lfc = os.path.join(result_path,'{analysis}','DEA_FILTERED_LFC.csv'),
    output:
        dea_lfc_heatmap = report(os.path.join(result_path,'{analysis}','plots','DEA_LFC_heatmap.png'),
                              caption="../report/lfc_heatmap.rst",
                              category="{}_dea_limma".format(config["project_name"]),
                              subcategory="{analysis}"),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/heatmap.yaml"
    log:
        os.path.join("logs","rules","lfc_heatmap_{analysis}.log"),
    params:
        partition=config.get("partition"),
    script:
        "../scripts/heatmap.R"