import argparse
import os
import pandas as pd
import pickle
import gzip
from sklearn.ensemble import RandomForestClassifier, VotingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GridSearchCV, train_test_split, cross_val_score
from sklearn.metrics import roc_auc_score

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

class MANAscore:
    def __init__(self):
        # Define internal file paths for training data
        self.file_paths = [
            os.path.join(BASE_DIR, '3gene/training/p2_known.csv'),
            os.path.join(BASE_DIR, '3gene/training/p2_known_RNA.csv'),
            os.path.join(BASE_DIR, '3gene/training/p11_known.csv'),
            os.path.join(BASE_DIR, '3gene/training/p11_known_RNA.csv'),
            os.path.join(BASE_DIR, '3gene/training/p15_known.csv'),
            os.path.join(BASE_DIR, '3gene/training/p15_known_RNA.csv')
        ]
        # Initialize attributes for models, classifiers, and data splits
        self.LM_models = []
        self.RF_models = []
        self.voting_i_classifier = None
        self.voting_ni_classifier = None
        # Initialize data split lists as instance attributes
        self.X_train_list = []
        self.X_test_list = []
        self.y_train_list = []
        self.y_test_list = []


    def load_and_split_data(self):
        """Load data and perform train-test split for each dataset internally."""
        self.X_train_list, self.X_test_list, self.y_train_list, self.y_test_list = [], [], [], []
        for path in self.file_paths:
            dat = pd.read_csv(path, index_col=0)
            X_train, X_test, y_train, y_test = train_test_split(
                dat[['CXCL13', 'ENTPD1', 'IL7R']],
                dat['Label'],
                test_size=0.2,
                random_state=42
            )
            self.X_train_list.append(X_train)
            self.X_test_list.append(X_test)
            self.y_train_list.append(y_train)
            self.y_test_list.append(y_test)

    def train_logistic_models(self):
        """Train logistic regression models for imputed and non-imputed classifiers."""
        LM_labels = ['LMi_p2', 'LMni_p2', 'LMi_p11', 'LMni_p11', 'LMi_p15', 'LMni_p15']
        for i, label in enumerate(LM_labels):
            model = LogisticRegression()
            model.fit(self.X_train_list[i], self.y_train_list[i])
            scores = cross_val_score(model, self.X_train_list[i], self.y_train_list[i], cv=5, scoring='roc_auc')
            auc_score = roc_auc_score(self.y_train_list[i], model.predict_proba(self.X_train_list[i])[:, 1])
            self.LM_models.append((label, model, scores.mean(), auc_score))

    def train_random_forest_models(self):
        """Train random forest models with GridSearchCV for imputed and non-imputed classifiers."""
        RF_labels = ['RFi_p2', 'RFni_p2', 'RFi_p11', 'RFni_p11', 'RFi_p15', 'RFni_p15']
        for i, label in enumerate(RF_labels):
            param_distributions = {'max_features': [1, 2, 3], 'n_estimators': [100, 500, 1000, 2000, 2500]}
            dt = RandomForestClassifier(random_state=42)
            grid = GridSearchCV(estimator=dt, param_grid=param_distributions, cv=10)
            grid.fit(self.X_train_list[i], self.y_train_list[i])
            best_model = grid.best_estimator_
            auc_score = roc_auc_score(self.y_train_list[i], best_model.predict_proba(self.X_train_list[i])[:, 1])
            self.RF_models.append((label, best_model, grid.best_params_, auc_score))

    def create_and_fit_voting_classifiers(self):
        """Create and fit voting classifiers for imputed and non-imputed models."""
        
        # Extract only label and model pairs for creating voting classifiers
        lm_dict = {label: model for label, model in self.LM_models}
        rf_dict = {label: model for label, model in self.RF_models}

        # Create and fit the voting classifier for imputed data
        self.voting_i_classifier = VotingClassifier(
            estimators=[
                ('LMi_p2', lm_dict['LMi_p2']),
                ('LMi_p11', lm_dict['LMi_p11']),
                ('LMi_p15', lm_dict['LMi_p15']),
                ('RFi_p2', rf_dict['RFi_p2']),
                ('RFi_p11', rf_dict['RFi_p11']),
                ('RFi_p15', rf_dict['RFi_p15'])
            ],
            voting='soft'
        )
        
        # Combine all imputed training data (indices 0, 2, 4) and labels
        X_train_i = pd.concat([self.X_train_list[idx] for idx in [0, 2, 4]])
        y_train_i = pd.concat([self.y_train_list[idx] for idx in [0, 2, 4]])
        
        # Fit the imputed classifier
        self.voting_i_classifier.fit(X_train_i, y_train_i)

        # Create and fit the voting classifier for non-imputed data
        self.voting_ni_classifier = VotingClassifier(
            estimators=[
                ('LMni_p2', lm_dict['LMni_p2']),
                ('LMni_p11', lm_dict['LMni_p11']),
                ('LMni_p15', lm_dict['LMni_p15']),
                ('RFni_p2', rf_dict['RFni_p2']),
                ('RFni_p11', rf_dict['RFni_p11']),
                ('RFni_p15', rf_dict['RFni_p15'])
            ],
            voting='soft'
        )
        
        # Combine all non-imputed training data (indices 1, 3, 5) and labels
        X_train_ni = pd.concat([self.X_train_list[idx] for idx in [1, 3, 5]])
        y_train_ni = pd.concat([self.y_train_list[idx] for idx in [1, 3, 5]])
        
        # Fit the non-imputed classifier
        self.voting_ni_classifier.fit(X_train_ni, y_train_ni)
    
    def fit_models(self):
        # Combine all imputed training data (indices 0, 2, 4) and labels
        X_train_i = pd.concat([self.X_train_list[idx] for idx in [0, 2, 4]])
        y_train_i = pd.concat([self.y_train_list[idx] for idx in [0, 2, 4]])

        # Fit the imputed classifier
        self.voting_i_classifier.fit(X_train_i, y_train_i)

        # Combine all non-imputed training data (indices 1, 3, 5) and labels
        X_train_ni = pd.concat([self.X_train_list[idx] for idx in [1, 3, 5]])
        y_train_ni = pd.concat([self.y_train_list[idx] for idx in [1, 3, 5]])

        # Fit the non-imputed classifier
        self.voting_ni_classifier.fit(X_train_ni, y_train_ni)
    
    def save_voting_models(self):
        """Save the fitted voting classifiers to compressed .pkl.gz files."""
        with gzip.open('voting_i_classifier.pkl.gz', 'wb') as f:
            pickle.dump(self.voting_i_classifier, f)
        with gzip.open('voting_ni_classifier.pkl.gz', 'wb') as f:
            pickle.dump(self.voting_ni_classifier, f)

    def load_voting_models(self):
        """Load the compressed .pkl.gz voting classifiers."""
        with gzip.open(os.path.join(BASE_DIR, 'models/voting_i_classifier.pkl.gz'), 'rb') as f:
            self.voting_i_classifier = pickle.load(f)
        with gzip.open(os.path.join(BASE_DIR, 'models/voting_ni_classifier.pkl.gz'), 'rb') as f:
            self.voting_ni_classifier = pickle.load(f)

    def evaluate_classifiers(self):
        """Evaluate both voting classifiers on combined imputed and non-imputed test data and return AUC scores."""
        if self.voting_i_classifier is None or self.voting_ni_classifier is None:
            raise ValueError("Voting classifiers must be created before evaluation. Call create_voting_classifiers() first.")
    
        # Check if test lists have enough data for expected indices
        required_indices_i = [0, 2, 4]
        required_indices_ni = [1, 3, 5]
        print(self.X_test_list)
        print(self.y_test_list)
        # Verify that required indices exist in X_test_list and y_test_list
        if not all(idx < len(self.X_test_list) for idx in required_indices_i + required_indices_ni):
            print(f"Error: One or more required test datasets are missing.")
            print(f"X_test_list length: {len(self.X_test_list)}, expected indices: {required_indices_i + required_indices_ni}")
            print(f"y_test_list length: {len(self.y_test_list)}, expected indices: {required_indices_i + required_indices_ni}")
            return None, None
    
        # Combine all imputed test data (indices 0, 2, 4) and labels
        X_test_i = pd.concat([self.X_test_list[idx] for idx in required_indices_i])
        y_test_i = pd.concat([self.y_test_list[idx] for idx in required_indices_i])
    
        # Combine all non-imputed test data (indices 1, 3, 5) and labels
        X_test_ni = pd.concat([self.X_test_list[idx] for idx in required_indices_ni])
        y_test_ni = pd.concat([self.y_test_list[idx] for idx in required_indices_ni])
    
        # Evaluate the combined imputed data with the imputed voting classifier
        try:
            auc_i = roc_auc_score(y_test_i, self.voting_i_classifier.predict_proba(X_test_i)[:, 1])
        except Exception as e:
            print(f"Error evaluating imputed classifier: {e}")
            auc_i = None

        # Evaluate the combined non-imputed data with the non-imputed voting classifier
        try:
            auc_ni = roc_auc_score(y_test_ni, self.voting_ni_classifier.predict_proba(X_test_ni)[:, 1])
        except Exception as e:
            print(f"Error evaluating non-imputed classifier: {e}")
            auc_ni = None
    
        print(f"Imputed Classifier AUC on Combined Data: {auc_i}")
        print(f"Non-Imputed Classifier AUC on Combined Data: {auc_ni}")
    
        return auc_i, auc_ni


    def predict_and_save_both(self, input_file_i, input_file_ni, output_file_i, output_file_ni):
        """Predict scores using both classifiers with two input files and save to separate CSV files."""

        # Load input data for imputation classifier prediction
        pdat_i = pd.read_csv(input_file_i, index_col=0)
        X_i = pdat_i[['CXCL13', 'ENTPD1', 'IL7R']]

        # Load input data for non-imputation classifier prediction
        pdat_ni = pd.read_csv(input_file_ni, index_col=0)
        X_ni = pdat_ni[['CXCL13', 'ENTPD1', 'IL7R']]

        # Check if voting classifiers are loaded
        if self.voting_i_classifier is None or self.voting_ni_classifier is None:
            raise ValueError("Both voting classifiers must be loaded. Call load_voting_models() first.")

        # Predictions using voting_i_classifier (imputed)
        prob_i = self.voting_i_classifier.predict_proba(X_i)[:, 1]
        d_i = pd.DataFrame({'barcode': pdat_i.index, 'score': prob_i})
        d_i.to_csv(output_file_i, index=False)
        print(f"Imputed classifier predictions saved to {output_file_i}")

        # Predictions using voting_ni_classifier (non-imputed)
        prob_ni = self.voting_ni_classifier.predict_proba(X_ni)[:, 1]
        d_ni = pd.DataFrame({'barcode': pdat_ni.index, 'score': prob_ni})
        d_ni.to_csv(output_file_ni, index=False)
        print(f"Non-imputed classifier predictions saved to {output_file_ni}")

def main():
    parser = argparse.ArgumentParser(description="MANAscore Tool")
    subparsers = parser.add_subparsers(dest="command", help="Available commands")
        
    train_parser = subparsers.add_parser("train", help="Train and save models")
    train_parser.add_argument("--save_models", action="store_true", help="Save trained models after training")
        
    load_parser = subparsers.add_parser("load", help="Load saved models")
    predict_parser = subparsers.add_parser("predict", help="Predict MANAscore using loaded models")
    predict_parser.add_argument("input_file_i", help="Path to input file for imputed classifier")
    predict_parser.add_argument("input_file_ni", help="Path to input file for non-imputed classifier")
    predict_parser.add_argument("output_file_i", help="Path to save output file for imputed classifier")
    predict_parser.add_argument("output_file_ni", help="Path to save output file for non-imputed classifier")
        
    args = parser.parse_args()
    manascore = MANAscore()
    if args.command == "train":
        manascore.load_and_split_data()
        manascore.train_logistic_models()
        manascore.train_random_forest_models()
        manascore.create_and_fit_voting_classifiers()
        
        
        if args.save_models:
            manascore.save_voting_models()
            print("Models saved successfully.")

    elif args.command == "load":
        print("Loading saved models...")
        manascore.load_voting_models()
        print("Models loaded successfully.")
    
    elif args.command == "predict":
        print(f"Predicting scores for {args.input_file_i} and {args.input_file_ni}...")
        if manascore.voting_i_classifier is None or manascore.voting_ni_classifier is None:
            manascore.load_voting_models()
        
        manascore.predict_and_save_both(args.input_file_i, args.input_file_ni, args.output_file_i, args.output_file_ni)
    else:
        parser.print_help()
if __name__ == "__main__":
    main()


