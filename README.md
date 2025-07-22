# BCH-ED
A repository of the code used to produce the results in our submission to PSB.

# Overview of files

NOTE: The user will have to specify directions for save_filepath in order to run these.

'preprocessing.r': cleans the original data, extracting demographic characteristics, chief complaints, and history of previous visits.

'produce-odds-ratios.r': uses the files produced in preprocessing.r to calculate odds ratios for admission and our secondary outcomes.

'make-forest-plots.r': produces forest plots for the paper.

'sample-visits-for-model-training.r': 'samples' each visit every 30 minutes and splits test/train data chronologically to train a predictive model. NOTE: this requires

'hyperparameter-tuning.py': uses 5-fold cross-validation to select hyperparameters on the training data.

'train-model-and-pull-shap-values.py': trains the model with these hyperparameters, makes predictions, and pulls Shapley values.

'make-beeswarm-plot.py': makes beeswarm plot from calculated SHAP values.

'plot-relative-shap-values.r': makes SHAP plots.
