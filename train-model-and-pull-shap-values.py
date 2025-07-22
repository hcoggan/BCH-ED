import numpy as np
import xgboost as xgb
import pandas as pd
from scipy import stats
from sklearn.model_selection import RandomizedSearchCV, KFold
from sklearn.metrics import roc_auc_score, average_precision_score, f1_score
import shap
import random

load_path = ""
save_path = ""

training_window_name = "four_years"
num_folds = 5 #5-fold cross-validation

train_data = pd.read_csv(load_path+ "training_window_"+ training_window_name+ "_training_data.csv", index_col=0)
test_data = pd.read_csv(load_path+ "training_window_"+ training_window_name+ "_test_data.csv", index_col=0)


#assign training weights- sum to 1 across a visit
train_data['timestamp_count'] = train_data.groupby('csn')['csn'].transform('count')
train_data['training_weight'] = 1 / train_data['timestamp_count']
train_data = train_data.drop(columns='timestamp_count')

test_data['timestamp_count'] = test_data.groupby('csn')['csn'].transform('count')
test_data['training_weight'] = 1 / test_data['timestamp_count']
test_data = test_data.drop(columns='timestamp_count')


column_labels = [name for name in train_data.columns.tolist() + test_data.columns.tolist() if name.lower().startswith("x")]

features_to_exclude = ["csn", "arrival_timestamp", "is_admitted", "is_discharged", "training_weight", "ed_complaint", "pediatric_comorbidity_score",
                             "mean_temperature", "max_temperature", "min_temperature", "mean_mean_arterial_pressure_device", "max_mean_arterial_pressure_device",
                             "min_mean_arterial_pressure_device"] + column_labels
Xs_train = train_data.drop(columns=features_to_exclude, errors="ignore") #ignore if cols don't exist
ys_train = train_data['is_admitted']
ids_train = train_data[['csn', 'timestamp']] #so we can identify this later
weights_train = train_data['training_weight'] #so all visits add up to 1 patient

Xs_test = test_data.drop(columns=features_to_exclude, errors="ignore") #ignore if cols don't exist
ys_test = test_data['is_admitted']
ids_test = test_data[['csn', 'timestamp']] #so we can identify this later
weights_test = test_data['training_weight'] #so all visits add up to 1 patient



best_params = {'n_estimators': 992,
              'learning_rate':  0.5974666765171989,
              'subsample': 0.3988175284874041,
              'max_depth': 9,
              'colsample_bytree': 0.7506276501255351,
              'min_child_weight': 3,
              'tree_method': ['gpu_hist'],
              'predictor': ['gpu_predictor']
             }



#parameters from CV
clf_xgb = xgb.XGBClassifier(objective = 'binary:logistic', eval_metric='logloss',
    use_label_encoder=False,
    random_state=2602,  gpu_id=0, params=best_params)

#train the model
clf_xgb.fit(Xs_train, ys_train, sample_weight=weights_train)

#make predictions on test data and save (this should get probability)
ys_pred = clf_xgb.predict_proba(Xs_test)[:, 1] #get the second column
ys_predicted_class = clf_xgb.predict(Xs_test) #get actual classifications

prediction_df = ids_test.copy()
prediction_df['ys_pred'] = ys_pred
prediction_df['ys_true'] = ys_test
prediction_df['weight'] = weights_test

prediction_df.to_csv(save_path+'time_based_predictions_on_test_data.csv', index=False)

#test these predictions
test_auroc = roc_auc_score(ys_test, ys_pred, sample_weight=weights_test)
test_precision_score = average_precision_score(ys_test, ys_pred, sample_weight=weights_test)
test_f1 = f1_score(ys_test, ys_predicted_class, sample_weight=weights_test)

print("AUROC on test data:", test_auroc)
print("Average precision score on test data:", test_precision_score)
print("F1 score on test data:", test_f1)

#now calculate Shapley values on 1000 predictions from test set
idx = random.sample(range(len(Xs_test)), 1000)

Xs_sampled = Xs_test.iloc[idx]
csns_sampled = ids_test.iloc[idx]

print("Columns used to train:")
print(Xs_sampled.columns.tolist())

explainer = shap.TreeExplainer(clf_xgb)
shap_values = explainer.shap_values(Xs_sampled)


feature_names = Xs_sampled.columns
shap_df = pd.DataFrame(shap_values, columns=feature_names)

csns_sampled = csns_sampled.reset_index(drop=True)
shap_df = shap_df.reset_index(drop=True)

#save these with csns as indices
shap_df = pd.concat([csns_sampled, shap_df], axis=1)
shap_df.to_csv(save_path+'time_based_shap_values.csv', index=False)
