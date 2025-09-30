# BCH-ED
A repository of the code used to produce the results in the paper

Coggan, H., Bischops, A., Chaudhari, P., Barak-Corren, Y., Fine, A. M., Reis, B. Y., & Aysola, J. (2025). Deciphering the influence of demographic factors on the treatment of pediatric patients in the emergency department. _Pacific Symposium on Biocomputing (PSB)_


# Overview of files

## Requirements

- R version 4.4.2 or higher
- `packages.txt` contains the R libraries that are needed. To install them, run `Rscript install-packages.r`.

## Usage

1. `bch-preprocess-data.r`: Preprocesses raw files, categorises and cleans variables, and attaches vital signs across stay.`
2. `bch-cohort-descriptors.r`: Produces the cohort descriptions used in Table 1 of the manuscript.
3. `bch-process-interventions.r`: Links interventions that occur during stay in the general terms required for the epidemiological analysis (whether or not a patient received any labs, any radiological tests, etc.)
4. `bch-get-interventions-for-prediction.r`: Links interventions that occur during stay in the *specific* terms required for the prediction study (e.g. whether and when a patient received any of the top 50 medication-route combinations, age-normalised lab results, etc), and produces training and test datasets.
5. `bch-produce-odds-ratios.r`: Runs the epidemiological analyses described in the paper.
6. `bch-make-forest-plots.r`: Produces forest plots to describe the epidemiological analyses.
7. `hyperparameter-tuning.py`: Conducts hyperparameter tuning to predict admission with 5-fold cross-validation.
8. `bch-get-shap.py`: Trains the model with chosen hyperparameters, and calculates Shapley values.
9. `make-shap-plot.ipynb`: Produces a beeswarm plot (Fig 5).
10. `bch-analyse-predictions.r`: Calculates AUROC and AUPRCs, as well as feature importance disparities.

# Acknowledgements

This research was supported by the National Library of Medicine of the National Institutes of Health (NIH) Award R01LM014300. 
The content is the responsibility of the authors and does not necessarily represent the official views of NIH. 
