import pandas as pd
import numpy as np
from sklearn.preprocessing import LabelEncoder, MinMaxScaler
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from sklearn.model_selection import train_test_split, StratifiedKFold, cross_val_score
from sklearn.pipeline import Pipeline
from sklearn.metrics import make_scorer, precision_score, accuracy_score, recall_score, f1_score, classification_report, confusion_matrix
from imblearn.over_sampling import BorderlineSMOTE
from xgboost import XGBClassifier


#load the data
data = pd.read_excel("Fig_2.xlsx", header=0)

x = data.iloc[:, 5:6].values  # Feature matrix
y = data.iloc[:, 0].values   # Target

print(y)

#label encoding
le = LabelEncoder()
y_label = le.fit_transform(y)

print(y_label)

#scaling the feature matrix
scaler_xlsx = MinMaxScaler(feature_range=(0, 1))  # LDA doesn't require uint8
X_scaled = scaler_xlsx.fit_transform(x)

#splitting the train and test set
X_train, X_test, y_train, y_test = train_test_split(X_scaled, y_label, test_size=0.4, stratify=y_label, random_state=42)
print(y_train)

#Apply Borderline SMOTE
min_class_count = min([np.sum(y_train == i) for i in np.unique(y_train)])

# Now, use this to set k_neighbors in SMOTE, considering there should be at least 2 samples 
# (SMOTE's requirement) to use as the nearest neighbors.
k_neighbors = max(min_class_count - 1, 1)  # Ensuring k_neighbors is at least 1

smote = BorderlineSMOTE(k_neighbors=k_neighbors)
X_train_smote, y_train_smote = smote.fit_resample(X_train, y_train)

#Pipeline for XGBoost and LDA
#pipeline = Pipeline(
#    [('lda', LDA()),
#     ('xgb', XGBClassifier(objective='multi:softprob', eval_metric='logloss', random_state=42, max_depth =3, min_child_weight=1, gamma=0, subsample=0.8, colsample_bytree=0.8, reg_lambda=1.0, reg_alpha = 0.5))]
#)

xgb = XGBClassifier(objective='multi:softprob', eval_metric='logloss', random_state=42, max_depth =3, min_child_weight=1, gamma=0, subsample=0.8, colsample_bytree=0.8, reg_lambda=1.0, reg_alpha = 0.5)

# Perform cross-validation
# Define scoring metrics
scoring = {
    'accuracy': make_scorer(accuracy_score),
    'precision_macro': make_scorer(precision_score, average='macro'),
    'recall_macro': make_scorer(recall_score, average='macro'),
    'f1_macro': make_scorer(f1_score, average='macro')
}

# Placeholder for true and predicted labels
true_labels = []
predicted_labels = []

# Set up StratifiedKFold
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

for train_idx, test_idx in cv.split(X_train_smote, y_train_smote):
    X_cv_train, X_cv_test = X_train_smote[train_idx], X_train_smote[test_idx]
    y_cv_train, y_cv_test = y_train_smote[train_idx], y_train_smote[test_idx]
    
    # Train the model on the cross-validation train set
    xgb.fit(X_cv_train, y_cv_train)
    
    # Predict on the cross-validation test set
    y_cv_pred = xgb.predict(X_cv_test)
    
    # Collect the true and predicted labels
    true_labels.extend(y_cv_test)
    predicted_labels.extend(y_cv_pred)

# Now, you can use classification_report on the collected labels
report = classification_report(true_labels, predicted_labels, target_names=le.classes_)
print(report)

cm = confusion_matrix(true_labels, predicted_labels)
print(cm)

#Evaluating model on testing set
# Train the model
xgb.fit(X_train_smote, y_train_smote)

# Make predictions
y_pred = xgb.predict(X_test)

# Calculate metrics
accuracy = accuracy_score(y_test, y_pred)
precision = precision_score(y_test, y_pred, average='macro')
recall = recall_score(y_test, y_pred, average='macro')
f1 = f1_score(y_test, y_pred, average='macro')

# Print the evaluation metrics
print(f"Accuracy: {accuracy:.4f}")
print(f"Precision (Macro): {precision:.4f}")
print(f"Recall (Macro): {recall:.4f}")
print(f"F1-Score (Macro): {f1:.4f}")

# Optionally, print the confusion matrix
conf_matrix = confusion_matrix(y_test, y_pred)
print("Confusion Matrix:")
print(conf_matrix)

 
