##### utility functions #####

def get_paths(wildcards):
    return [annot.loc[wildcards.analysis,'data'], annot.loc[wildcards.analysis,'metadata']]

def get_feature_list_path(wildcards):
    return config["feature_lists"][wildcards.feature_list]

def get_feature_lists(wildcards):
    checkpoint_output = checkpoints.aggregate.get(**wildcards).output.feature_lists
    return expand(os.path.join(checkpoint_output, "{i}"),
                  i=glob_wildcards(os.path.join(checkpoint_output, "{i}")).i)