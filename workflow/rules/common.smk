##### utility functions #####

def get_paths(wildcards):
    return [annot.loc[wildcards.analysis,'data'], annot.loc[wildcards.analysis,'metadata']]
