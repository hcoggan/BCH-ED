import numpy as np
import xgboost as xgb
import pandas as pd
from scipy import stats
from sklearn.model_selection import RandomizedSearchCV, KFold
from sklearn.metrics import roc_auc_score
import argparse

load_path= ""
save_path = ""

num_folds = 5 #5-fold cross-validation


#Load data
train_data = pd.read_csv(load_path+ "training-data.csv", index_col=0)
              

#assign training weights- sum to 1 across a visit
train_data['timestamp_count'] = train_data.groupby('csn')['csn'].transform('count')
train_data['training_weight'] = 1 / train_data['timestamp_count']
train_data = train_data.drop(columns='timestamp_count')

#Exclude counting labels, which start with X or V.
column_labels = [name for name in train_data.columns.tolist() if (name.lower().startswith("x") or name.lower().startswith("v"))]


features_to_exclude = ["csn", "arrival_timestamp", "is_admitted", "training_weight"] + column_labels
Xs_train = train_data.drop(columns=features_to_exclude, errors="ignore") #ignore if cols don't exist
ys_train = train_data['is_admitted']
ids_train = train_data[['csn', 'timestamp']] #so we can identify this later
weights_train = train_data['training_weight'] #so all visits add up to 1 patient



#parameters from SO 43927725
clf_xgb = xgb.XGBClassifier(objective = 'binary:logistic', eval_metric='logloss',
    use_label_encoder=False,
    random_state=2602, tree_method='gpu_hist', gpu_id=0)

param_dist = {'n_estimators': stats.randint(150, 1000),
              'learning_rate': stats.uniform(0.01, 0.59),
              'subsample': stats.uniform(0.3, 0.6), #searches 0.3 to 0.9
              'max_depth': [3, 4, 5, 6, 7, 8, 9],
              'colsample_bytree': stats.uniform(0.5, 0.4), #searches 0.5 to 0.9
              'min_child_weight': [1, 2, 3, 4],
              'tree_method': ['gpu_hist'],
              'predictor': ['gpu_predictor']
             }


kfold_5 = KFold(n_splits=num_folds, shuffle = True, random_state=2602)


def cv(num_iterations):

    clf = RandomizedSearchCV(clf_xgb, 
                            param_distributions = param_dist,
                            cv = kfold_5,  
                            n_iter = num_iterations, 
                            scoring = 'roc_auc',
                            n_jobs=-1,
                            verbose=2)


    clf.fit(Xs_train, ys_train, sample_weight=weights_train)

    # Output best parameters and score
    print("Best parameters found:", clf.best_params_)
    print("Best cross-validated ROC AUC:", clf.best_score_)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Cross-validation of model trained on BCH data')
    parser.add_argument('num_iter', 
                       type=int,  
                       help='Number of iterations to use when hyperparameter tuning')
    args = parser.parse_args()
    num_iter = args.num_iter
    cv(num_iter)