# Configuration

You need one configuration file and one annotation file to run the complete workflow. You can use the provided example as starting point. Always use absolute paths. If in doubt read the comments in the config and/or try the default values.

- project configuration (config/config.yaml): different for every project/dataset and configures the analyses to be performed
- sample annotation (sample_annotation): CSV file consisting of 10 mandatory columns
    -  name: name of the dataset/analysis (tip: keep it short, but descriptive, distinctive and unique)
    -  data: absolute path to the input data as CSV file as feature by sample table (eg count matrix) that has already been quality controlled (eg bad samples removed) and filtered (eg only expressed genes)
    -  metadata: absolute path to the metadata for the required analysis (ie every variable in the formula needs to have a corresponding column)
    -  formula: ---
    -  block_var: ----
    -  comparisons: ---
    -  calcNormFactors_method: ---
    -  voom: ---
    -  eBayes: ---
    -  limma_trend: ---

