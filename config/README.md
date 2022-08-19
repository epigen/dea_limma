# Configuration

You need one configuration file and one annotation file to run the complete workflow. Additionally, you can provide feature annotation file in the project configuration (eg for plotting gene symbols instead of ensembl terms). You can use the provided example as starting point. Always use absolute paths. If in doubt read the comments in the config and/or try the default values.

- project configuration (config/config.yaml): configures the analyses to be performed and is different for every project/dataset.
- annotation (annotation): CSV file consisting of 10 mandatory columns
    -  name: name of the dataset/analysis (tip: keep it short, but descriptive, distinctive and unique)
    -  data: Absolute path to the input data as CSV file as feature by sample table (eg RNA count matrix) that has already been quality controlled (eg bad samples removed) and filtered for relevant features (eg only expressed genes). The first column has to contain the features and the first row the sample-names.
    -  metadata: Absolute path to the metadata as CSV file for the required analysis (ie every variable in the formula and opyional blocking variable needs to have a corresponding column). The first column has to be the sample name. The metadata file has to be R compatible (eg column names should not start with a number or contain colons).
    -  formula: A string that will be converted to a formula in R (eg ~ treatment + batch).
    -  block_var: A variable (present in the metadata as column) for the blocking feature (see README > Features).
    -  comparisons: Variable names contained in the formula (and metadata) which coefficient's you are interested in, separated by '|' (eg treatment|batch). Results of all derived groups (eg treatmentLP) containing one of the comparisons will be returned.
    -  calcNormFactors_method: Flag to indicate if [edgeR:calcNormFactors](https://www.rdocumentation.org/packages/edgeR/versions/3.14.0/topics/calcNormFactors) function be used specifing the parameter "method" (eg none or TMM) or should be skipped (0).
    -  voom: Flag to indicate if voom function should be used (1) or not (0).
    -  eBayes: Flag to indicate if eBayes function should be used (1) or not (0).
    -  limma_trend: Flag to indicate if limma-trend should be used (1) (ie sets [limma::eBayes](https://www.rdocumentation.org/packages/limma/versions/3.28.14/topics/ebayes) parameter trend=TRUE), or not (0). Please make sure to activate the required eBayes function (=1) and deactivate voom (=0) if you use limma-trend. Using voom and limma-trend makes no sense, but is not forbiden by the workflow.

