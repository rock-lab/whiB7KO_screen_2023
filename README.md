# whiB7KO screen (2023)
Repository containing the code and methods for Poulton et al. 2023


Below is a short legend describing the main files, followed by example running the code:

#### vulnerability_modeling_tools.R
R code with functions useful for running the vulnerability analysis.

#### vulnerability_plotting_tools.R
R code with functions useful for plotting vulnerability results

#### gene_vulnerability_analysis.R
R code for running vulnerability analysis in paralell.

#### model_negbinom_2line_logistic_w_guide_lambda_betamin_zero.stan
Stan code for the vulnerability model.

#### gene_vulnerability_plots.R
R code for creating plots of vulnerability results in parallel.



# Example

### Analysis
Users can run the vulnerability analyis performed in the paper using the following script:

`gene_vulnerability_analysis.R`

The following parameters must be passed:

 - **data**:  path to dataframe containing the counts data
 - **strain**: strain of the organism (for tracking/labeling purposes) 
 - **label**: short label (for tracking/labeling purposes)
 - **output**: output directory for the samples and model data for each gene
 - **cores**: number of cores to run in parallel


Example:

`Rscript gene_vulnerability_analysis.R --data data/data_counts_atc_Smeg_whiB7KO_08_10_2023.txt --strain Smeg --label Smeg_whiB7KO --output results/ --cores 40`


### Plotting

Similarly, users can create plots of the results by running the following script:

`gene_vulnerability_plots.R`

Users can provide four parameters:

 - **data**:  path to dataframe containing the counts data
 - **strain**: strain of the organism (for tracking/labeling purposes) 
 - **exp**: short label for the experiment (for tracking/labeling purposes)
 - **date**: date where the run, experiments, or processing were done (for tracking/labeling purposes)
 - **plot_dir**: output directory for the plots
 - **dir**: directory where the samples and model data were outputed by the sampling process
 - **cores**: number of cores to run in parallel

Example:

`Rscript gene_vulnerability_plots.R --data data/data_counts_atc_Smeg_whiB7KO_08_10_2023.txt --strain Smeg --exp whiB7KO --date 08_10_2023 --plot_dir plots/ --dir results/ --cores 10`





