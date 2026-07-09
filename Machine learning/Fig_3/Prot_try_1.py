import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import LabelEncoder, MinMaxScaler
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split, StratifiedKFold, cross_val_score
from sklearn.pipeline import Pipeline
from sklearn.metrics import make_scorer, precision_score, accuracy_score, recall_score, f1_score, classification_report, confusion_matrix, roc_auc_score, roc_curve
from imblearn.over_sampling import BorderlineSMOTE
from xgboost import XGBClassifier
from sklearn.metrics import roc_curve, roc_auc_score
from sklearn.preprocessing import label_binarize
from sklearn.multiclass import OneVsRestClassifier
import seaborn as sns


#load the data
data = pd.read_excel("Prot.xlsx", header=0)

x = data.iloc[:, 6:11].values  # Feature matrix
y = data.iloc[:, 0].values   # Target

# For the example, using column indices 5:6 as per your setup
feature_columns = data.columns[6:11]

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

# Initialize an empty dictionary to store AUC scores for each feature
feature_auc_scores = {}

# Initialize OneVsRestClassifier with XGBClassifier
ovr_clf = OneVsRestClassifier(XGBClassifier(use_label_encoder=False, eval_metric='mlogloss', random_state=42, max_depth=3, min_child_weight=1, gamma=0, subsample=0.8, colsample_bytree=0.8, reg_lambda=1.0, reg_alpha=0.5))
#ovr_clf = OneVsRestClassifier(LogisticRegression())

# Train the classifier with SMOTE-applied training data
ovr_clf.fit(X_train_smote, y_train_smote)

# Make predictions with the test set
y_pred = ovr_clf.predict(X_test)  # These are class labels, not probabilities

# Classification report
print(classification_report(y_test, y_pred, target_names=le.classes_))

# For ROC curve and AUC, iterate over each class
for i, class_name in enumerate(le.classes_):
    # Compute ROC curve and AUC for each class
    y_test_bin = label_binarize(y_test, classes=range(len(le.classes_)))[:, i]  # Binarize labels for the current class
    y_score = ovr_clf.predict_proba(X_test)[:, i]  # Get probabilities for the current class
    
    fpr, tpr, _ = roc_curve(y_test_bin, y_score)
    roc_auc = roc_auc_score(y_test_bin, y_score)

    # Plotting
    plt.figure()
    plt.plot(fpr, tpr, lw=2, label=f'Class {class_name} (AUC = {roc_auc:.2f})')
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(f'Receiver Operating Characteristic for Class {class_name}')
    plt.legend(loc="lower right")
    plt.show()


######## for features individually
for j, feature_name in enumerate(feature_columns):
    # Selecting a single feature for the model
    X_train_feature = X_train_smote[:, j:j+1]
    X_test_feature = X_test[:, j:j+1]
    
    # Initialize OneVsRest classifier with XGBClassifier for the current feature
    ovr_clf_1 = OneVsRestClassifier(XGBClassifier(use_label_encoder=False, eval_metric='mlogloss', random_state=42, max_depth=3, min_child_weight=1, gamma=0, subsample=0.8, colsample_bytree=0.8, reg_lambda=1.0, reg_alpha=0.5))
    #ovr_clf_1 = OneVsRestClassifier(LogisticRegression())
    
    ovr_clf_1.fit(X_train_feature, y_train_smote)
    
    # Predicting class labels for the test set
    y_pred = ovr_clf_1.predict(X_test_feature)
    
    # Predicting probabilities for ROC AUC
    y_proba = ovr_clf_1.predict_proba(X_test_feature)
    
    # Calculate and print the classification report for the current feature
    print(f"Classification report for {feature_name}:")
    print(classification_report(y_test, y_pred, target_names=le.classes_))


df = pd.read_excel("F1_prot.xlsx")

df.set_index('Feature_Name', inplace = True) 

# Assuming 'df' is your DataFrame prepared as above
plt.figure(figsize = (20,10))
sns.heatmap(df, annot=True, cmap='coolwarm', fmt=".2f")
plt.title('Predictive Classification')
plt.ylabel('Feature')
plt.xlabel('Class')
plt.savefig('Biomarkers_F1.jpeg')