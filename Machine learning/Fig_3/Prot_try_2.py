import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import LabelEncoder, MinMaxScaler
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split, StratifiedKFold, cross_val_score, cross_validate
from sklearn.pipeline import Pipeline
from sklearn.metrics import make_scorer, precision_score, accuracy_score, recall_score, f1_score, classification_report, confusion_matrix, roc_auc_score, roc_curve
from imblearn.over_sampling import BorderlineSMOTE
from xgboost import XGBClassifier
from sklearn.metrics import roc_curve, roc_auc_score
from sklearn.preprocessing import label_binarize
from sklearn.multiclass import OneVsRestClassifier

#load the data
data = pd.read_excel("Prot.xlsx", header=0)

x = data.iloc[:, 2:].values  # Feature matrix
y = data.iloc[:, 0].values   # Target

# For the example, using column indices 5:6 as per your setup
feature_columns = data.columns[2:]

#label encoding
le = LabelEncoder()
y_label = le.fit_transform(y)

# Now, determine the number of unique classes
n_classes = len(le.classes_)

# Now, determine the number of unique classes
n_classes = len(le.classes_)

#scaling the feature matrix
scaler_xlsx = MinMaxScaler(feature_range=(0, 1))  # LDA doesn't require uint8
X_scaled = scaler_xlsx.fit_transform(x)

#splitting the train and test set
X_train, X_test, y_train, y_test = train_test_split(X_scaled, y_label, test_size=0.4, stratify=y_label, random_state=42)

# Correctly calculate the minimum class count based on y_train (not binarized)
min_class_count = min(np.bincount(y_train))  # np.bincount efficiently counts occurrences of each value in an array of non-negative ints.

# Adjust k_neighbors based on the smallest class size for SMOTE
k_neighbors = max(min_class_count - 1, 1)

smote = BorderlineSMOTE(k_neighbors=k_neighbors)
X_train_smote, y_train_smote = smote.fit_resample(X_train, y_train)

#define pipeline
xgb = XGBClassifier(objective='multi:softprob', eval_metric='mlogloss', random_state=42, max_depth = 3, min_child_weight=1, gamma=0, subsample=0.8, colsample_bytree=0.8, reg_lambda=1.0, reg_alpha = 0.5)


cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

scoring = {
    'accuracy': make_scorer(accuracy_score),
    'precision_macro': make_scorer(precision_score, average='macro'),
    'recall_macro': make_scorer(recall_score, average='macro'),
    'f1_macro': make_scorer(f1_score, average='macro')
}

for j, feature_name in enumerate(feature_columns):
    # Selecting a single feature for the model
    X_train_feature = X_train_smote[:, j:j+1]
    X_test_feature = X_test[:, j:j+1]
    
    # Initialize OneVsRest classifier with XGBClassifier for the current feature
    
    xgb.fit(X_train_feature, y_train_smote)
    
    # Predicting class labels for the test set
    y_pred = xgb.predict(X_test_feature)
    
    # Predicting probabilities for ROC AUC
    y_proba = xgb.predict_proba(X_test_feature)
    
    # Calculate and print the classification report for the current feature
    print(f"Classification report for {feature_name}:")
    print(classification_report(y_test, y_pred, target_names=le.classes_))

    # Perform cross-validation
    cv_results = cross_validate(xgb, X_train_smote, y_train_smote, cv=cv, scoring=scoring)
    print(cv_results)

######## for features individually
for j, feature_name in enumerate(feature_columns):
    # Selecting a single feature for the model
    X_train_feature = X_train_smote[:, j:j+1]
    X_test_feature = X_test[:, j:j+1]
    
    # Initialize OneVsRest classifier with XGBClassifier for the current feature
    clf_1 = XGBClassifier(objective='multi:softprob', eval_metric='mlogloss', random_state=42, max_depth = 3, min_child_weight=1, gamma=0, subsample=0.8, colsample_bytree=0.8, reg_lambda=1.0, reg_alpha = 0.5)
    #ovr_clf_1 = OneVsRestClassifier(LogisticRegression())
    
    clf_1.fit(X_train_feature, y_train_smote)
    
    # Predicting class labels for the test set
    y_pred = clf_1.predict(X_test_feature)
    
    # Predicting probabilities for ROC AUC
    y_proba = clf_1.predict_proba(X_test_feature)
    
    # Calculate and print the classification report for the current feature
    print(f"Classification report for {feature_name}:")
    print(classification_report(y_test, y_pred, target_names=le.classes_))