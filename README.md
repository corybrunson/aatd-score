
# Risk scores for severe AATD

<!-- badges: start -->
<!-- badges: end -->

This repo contains code used in [the study by Riley, Brunson, Eydgahi, Brantly, and Lascano](https://openres.ersjournals.com/content/9/5/00302-2023).[^1]

## Objective

The aim of the study was to design an interpretable clinical decision support tool to inform decision to screen for abnormal AAT genotype.
The modeling goal was to predict abnormal genotype from history of hepatic disease, respiratory disease, and (optionally) smoking.

## Data source

This project relied on data collected and provided by the National Detection Program based at the University of Florida and GeneAidyx.

## Software dependencies

The project used [Tidyverse](https://www.tidyverse.org/) for general data science and [Tidymodels](https://www.tidymodels.org/) for machine learning.
Additional risk models were built using [FasterRisk](https://github.com/jiachangliu/FasterRisk).

Most scripts source project-wide settings from 'code/settings.r'.
Due to FasterRisk being implemented in Python and not accessible through Tidymodels, separate scripts were prepared to build and analyze these models.

## Model selection

The first round of analysis (the scripts '1-*.r') compared the performance of several families of predictive models on a subset of the data, using a simple train--test partition without hyperparameter optimization.
The goal was to inform decisions about a cross-validated modeling effort using the complete data set.
Specifically, we used this preliminary analysis to
(1) select predictors to consider,
(2) decide what interpretable model families to use, and
(3) how much predictive performance might be weakened by using interpretable rather than non-interpretable models.

## Model optimization

The second round of analysis (the scripts '2-*.r') conducted a stratified 6-fold cross-validated comparison of models selected based on the results of the first round:
(1) as predictors, either lung and liver history alone or together with either gender or smoking history; and
(2) as model families, logistic regression, random forest, nearest neighbor, and FasterRisk, together with a stand-in for the current guidance to screen patients with COPD.
This round used a majority training set of the data to determine optimal hyperparameter settings, with the remaining testing set held out for evaluation in the third round.

## Model evaluation

The third and final round of analysis (the scripts '3-*.r') fitted the optimized model families and specifications to the full training set from the second round and evaluated them on the testing set.

## Contact

Please email Cory Brunson (@corybrunson) with any question about the analysis and code.

[^1]: [Riley EL, Brunson JC, Eydgahi S, Brantly ML, Lascano JE. Development of a risk score to increase detection of severe alpha-1 antitrypsin deficiency. _ERJ Open Res_. 2023 Sep 18;9(5):00302-2023. doi: 10.1183/23120541.00302-2023. PMID: 37727673; PMCID: PMC10505949.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10505949/)
