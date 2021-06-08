from __future__ import absolute_import, division, print_function, unicode_literals
import pandas as pd
import numpy as np
import re
import pathlib
import seaborn as sns
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
from sklearn.metrics import classification_report, confusion_matrix, accuracy_score
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.model_selection import GridSearchCV
from sklearn.ensemble import RandomForestClassifier

def tuneHyperParam(cell_type, X_train, y_train):
    # define models and parameters
    model = RandomForestClassifier()
    n_estimators = [10, 50, 100, 500, 1000]
    max_features = ['sqrt', 'log2']
    # define grid search
    grid = dict(n_estimators=n_estimators,max_features=max_features)
    cv = RepeatedStratifiedKFold(n_splits=10, n_repeats=5, random_state=1)
    grid_search = GridSearchCV(estimator=model, param_grid=grid, n_jobs=-1, cv=cv, scoring='accuracy',error_score=0)
    grid_result = grid_search.fit(X_train, y_train)

    # Save Grid Search Results for Plotting later
    print("Best: %f using %s" % (grid_result.best_score_, grid_result.best_params_))
    
    X = pd.concat([pd.DataFrame(grid_result.cv_results_["params"]),
        pd.DataFrame(grid_result.cv_results_["mean_test_score"], columns=["Accuracy"]),
        pd.DataFrame(grid_result.cv_results_["std_test_score"], columns=["Std"])],axis=1)
    
    X.to_csv('./' + cell_type + '_tuneHyperParam.csv', sep=',',index=False)    
    
    return grid_result.best_params_['n_estimators'], grid_result.best_params_['max_features']

def rfc(counts_subset, cell_type, X_train, X_test, y_train, y_test, optimal_n_estimators, optimal_max_features):

    clf = RandomForestClassifier(n_estimators = optimal_n_estimators, max_features = optimal_max_features)

    # Train the model using the training sets y_pred=clf.predict(X_test)
    clf.fit(X_train,y_train)
    
    # Predict
    y_pred = clf.predict(X_test)
    y = pd.DataFrame()
    y["y_test"] = y_test
    y["y_pred"] = y_pred
    y.to_csv('./' + cell_type + '_predictions.csv', sep=',',index=True)
    print("Accuracy" + str(accuracy_score(y_test, y_pred)))
    
    # Get the feature importance scores
    feature_imp = pd.Series(clf.feature_importances_, index = list(counts_subset)[0:-1]).sort_values(ascending=False)
    feature_imp.to_csv('./' + cell_type + '_feature_imp.csv', sep=',',index=True)

    return

def runCellType(counts, df, cell_type):
    # Subset the counts matrix to get current cell type
    meta_data_subset = df[df['predicted.celltype.l2'] == cell_type]
    cells = list(meta_data_subset.index)
    counts_subset = counts.loc[cells,:]

    # Split the data into train and test datasets
    X = counts_subset[list(counts_subset)[0:-1]]  # Features
    y = counts_subset['condition']  # Labels

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3)

    # Get the optimal hyperparameters using the grid search method
    print("Optimizing Hyper-parameters...")
    optimal_n_estimators, optimal_max_features = tuneHyperParam(cell_type, X_train, y_train)

    # Run the random forest classification
    print("Running Random Forest Classification...")
    rfc(counts_subset, cell_type, X_train, X_test, y_train, y_test, optimal_n_estimators, optimal_max_features)

    return

if __name__ == '__main__':
    # Load in the Data
    print("Loading Data...")
    df = pd.read_csv('./alive_dead_day0_meta.csv', index_col=0)
    df = df.sort_index()
    sct_normalized_matrix = './alive_dead_day0_SCT_normalized.txt'
    counts = pd.read_csv(sct_normalized_matrix, sep=',', index_col=0).transpose()
    counts = counts.sort_index()
    counts["condition"] = np.array([i == "Alive_D0" for i in list(df["condition"])]).astype('int')
    cell_types = list(df['predicted.celltype.l2'].value_counts()[df['predicted.celltype.l2'].value_counts()>100].index)
    cell_types = list(reversed(cell_types))
    for cell_type in cell_types:
        print("Processing " + cell_type)
        runCellType(counts, df, cell_type)









