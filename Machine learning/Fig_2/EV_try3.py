import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder, MinMaxScaler
from sklearn.metrics import classification_report, precision_score, make_scorer, recall_score, f1_score, accuracy_score
from sklearn.pipeline import make_pipeline
from imblearn.over_sampling import BorderlineSMOTE
import numpy as np

# Models to compare
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from xgboost import XGBClassifier
from sklearn.model_selection import cross_val_score, StratifiedKFold



# Load data
data = pd.read_excel("Fig_2.xlsx", header=0)

# Separate features and target
X = data.iloc[:, 1:].values  # Adjust columns as needed
y = data.iloc[:, 0].values

# Encode target labels
le = LabelEncoder()
y_encoded = le.fit_transform(y)

# Scale features
scaler = MinMaxScaler()
X_scaled = scaler.fit_transform(X)

# Split data into training and test sets
X_train, X_test, y_train, y_test = train_test_split(X_scaled, y_encoded, test_size=0.4, stratify=y_encoded, random_state=42)

# Correctly calculate the minimum class count based on y_train (not binarized)
min_class_count = min(np.bincount(y_train))  # np.bincount efficiently counts occurrences of each value in an array of non-negative ints.

# Adjust k_neighbors based on the smallest class size for SMOTE
k_neighbors = max(min_class_count - 1, 1)

smote = BorderlineSMOTE(k_neighbors=k_neighbors)
X_train_smote, y_train_smote = smote.fit_resample(X_train, y_train)

models = {
    "Logistic Regression": LogisticRegression(max_iter=1000),
    "Random Forest": RandomForestClassifier(),
    "XGBoost": XGBClassifier(use_label_encoder=False, eval_metric='mlogloss')
}

# Evaluate each model
for name, model in models.items():
    print(f"Model: {name}")
    for i, feature_name in enumerate(data.columns[1:]): 
        X_train_feature = X_train_smote[:, i].reshape(-1, 1)
        X_test_feature = X_test[:, i].reshape(-1, 1)
        
        # Fit model
        model.fit(X_train_feature, y_train_smote)
        
        # Predict and evaluate
        y_pred = model.predict(X_test_feature)
        print(f"Feature: {feature_name}")
        print(classification_report(y_test, y_pred, target_names=le.classes_))


# Example with Logistic Regression for all features
model = XGBClassifier(use_label_encoder=False, eval_metric='mlogloss')
cv = StratifiedKFold(n_splits=5)

# Create a scorer object
accuracy_macro_scorer = make_scorer(accuracy_score)
precision_macro_scorer = make_scorer(precision_score, average='macro')
recall_macro_scorer = make_scorer(recall_score, average='macro')
f1_macro_scorer = make_scorer(f1_score, average='macro')

accuracy = {}
precision = {}
recall = {}
f1 = {}

for i, feature_name in enumerate(data.columns[1:]):  # Adjust based on your dataset's feature columns
    # Select the current feature
    X_train_feature = X_train_smote[:, i].reshape(-1, 1)
    
     # Perform cross-validation for the current feature
    accuracy_scores = cross_val_score(model, X_train_feature, y_train_smote, cv=cv, scoring= accuracy_macro_scorer)
    # Store the average score for the feature
    accuracy[feature_name] = np.mean(accuracy_scores)

    # Perform cross-validation for the current feature
    precision_scores = cross_val_score(model, X_train_feature, y_train_smote, cv=cv, scoring=precision_macro_scorer)
    # Store the average score for the feature
    precision[feature_name] = np.mean(precision_scores)

    # Perform cross-validation for the current feature
    recall_scores = cross_val_score(model, X_train_feature, y_train_smote, cv=cv, scoring=recall_macro_scorer)
    # Store the average score for the feature
    recall[feature_name] = np.mean(recall_scores)

    # Perform cross-validation for the current feature
    f1_scores = cross_val_score(model, X_train_feature, y_train_smote, cv=cv, scoring=f1_macro_scorer)
    # Store the average score for the feature
    f1[feature_name] = np.mean(f1_scores)

# Print the precision for each feature
for feature, score in accuracy.items():
    print(f"{feature}: Accuracy = {score:.3f}")

# Print the precision for each feature
for feature, score in precision.items():
    print(f"{feature}: Precision (Macro) = {score:.3f}")

# Print the recall for each feature
for feature, score in recall.items():
    print(f"{feature}: Recall (Macro) = {score:.3f}")

# Print the f1 for each feature
for feature, score in f1.items():
    print(f"{feature}: F1 (Macro) = {score:.3f}")

