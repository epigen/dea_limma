##### utility functions #####

def get_paths(wildcards):
    return [annot.loc[wildcards.analysis,'data'], annot.loc[wildcards.analysis,'metadata']]

def get_metadata_path(wildcards):
    return annot.loc[wildcards.analysis,'metadata']

def get_feature_list_path(wildcards):
    return config["feature_lists"][wildcards.feature_list]

def get_input_ova_stats_plot(wildcards):
    return [
        os.path.join(result_path, analysis,'results.csv')
        for analysis in std_analyses_to_ova_dict[wildcards.analysis]
    ]