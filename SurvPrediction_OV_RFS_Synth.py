import os
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error

from sklearn import set_config
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import OrdinalEncoder
from sksurv.util import Surv
from sklearn.model_selection import GridSearchCV

from sksurv.datasets import load_gbsg2
from sksurv.preprocessing import OneHotEncoder
from sksurv.ensemble import RandomSurvivalForest
import winsound

from sdv.single_table import CTGANSynthesizer
from sdv.metadata import SingleTableMetadata

import numpy as np
from skopt import BayesSearchCV
from skopt import gp_minimize
from sklearn.model_selection import cross_val_score
import joblib

# Setting up the Paths ############################################################################
directory = os.getcwd()

generate_synthetic_data = False
###################################################################################################

# Load the dataset
data = pd.read_csv(directory+"\\data\\surv_data.csv", sep = ",", decimal='.')
print(data.head())
print("Data is Loaded")

# Split the data for generating synthetic data and testing
data['group'] = pd.cut(data['OS_MONTHS'], bins=range(0, int(data['OS_MONTHS'].max())+5, 5), right=False)
#print(data.head())
# Filter the DataFrame to keep only the rows with groups that have at least 2 members
data = data[data['group'].isin(data['group'].value_counts()[data['group'].value_counts() >= 2].index)]
#print(data['group'].value_counts().sort_index())

# Finally make the split
data_train, data_test = train_test_split(data, test_size=0.15, random_state=42, stratify=data['group'])
data_train.drop(["group"], axis=1, inplace=True)
data_test.drop(["group"], axis=1, inplace=True)
print("Data is Split")

if generate_synthetic_data == True:
    # Create metadata from the data for synthetic generation
    metadata = SingleTableMetadata()
    metadata.detect_from_dataframe(data=data_train)
    print("Obtained Metadata")
    
    data_train_1 = data_train[data_train['SURV_STATUS'] == 1]
    data_train_0 = data_train[data_train['SURV_STATUS'] == 0]
    
    # Define Synthesizer
    nepochs = 150
    synthesizer_1 = CTGANSynthesizer(metadata, enforce_rounding=True, enforce_min_max_values = True, verbose = True, epochs = nepochs)
    synthesizer_0 = CTGANSynthesizer(metadata, enforce_rounding=True, enforce_min_max_values = True, verbose = True, epochs = nepochs)
    # Fit Synthesizer
    synthesizer_1.fit(data_train_1)
    synthesizer_0.fit(data_train_0)
    # Obtain Synthetic data
    nrows = 250
    synthetic_data_1 = synthesizer_1.sample(num_rows=nrows)
    synthetic_data_0 = synthesizer_0.sample(num_rows=nrows)
    
    synthetic_data = pd.concat([synthetic_data_1, synthetic_data_0], ignore_index=True)
    synthetic_data = synthetic_data.sample(frac=1, random_state=42)
    
    print(synthetic_data_1.size)
    print(synthetic_data_0.size)
    print(synthetic_data.size)
    
    synthetic_data.to_csv(directory+'\\data\\synthetic_data.csv', index=False)
else:
    synthetic_data = pd.read_csv(directory+"\\data\\synthetic_data.csv", sep = ",", decimal='.')

# split synthetic and testing data into X and y 
synthetic_data_X = synthetic_data.drop(["OS_MONTHS", "SURV_STATUS"], axis=1).copy()
synthetic_data_y = Surv.from_arrays(synthetic_data['SURV_STATUS'], synthetic_data['OS_MONTHS'])
data_test_X = data_test.drop(["OS_MONTHS", "SURV_STATUS"], axis=1).copy()
data_test_y = Surv.from_arrays(data_test['SURV_STATUS'], data_test['OS_MONTHS'])

# Define the hyperparameter search space
param_space = {
    'n_estimators': (500, 1000),
    'min_samples_split': (2, 50),
    'min_samples_leaf': (2, 30),
}

# Define the Bayesian optimizer
optimizer = BayesSearchCV(
    RandomSurvivalForest(n_jobs=-1, random_state=42),
    param_space,
    n_iter=50,  # Number of iterations for the Bayesian optimization
    n_points=1,  # Number of random points to probe the objective function before fitting a Gaussian Process
    cv=3,  # Cross-validation folds
    n_jobs=-1,
    random_state=42
)

# Fit the Bayesian optimizer to your data
optimizer.fit(synthetic_data_X, synthetic_data_y)

# Get the best parameters
best_params = optimizer.best_params_
#best_params = {'n_estimators':2500, "min_samples_split":2, "min_samples_leaf":2}
print("Best Parameters:", best_params)

# Use the best parameters to train your final model
rsf = RandomSurvivalForest(
    n_estimators=best_params['n_estimators'],
    min_samples_split=best_params['min_samples_split'],
    min_samples_leaf=best_params['min_samples_leaf'],
    n_jobs=-1,
    random_state=42
)

rsf.fit(synthetic_data_X, synthetic_data_y)

# Predict
y_pred = rsf.predict(data_test_X)
#y_pred = y_pred-50
y_actual = data_test['OS_MONTHS']

mse = mean_squared_error(y_actual, y_pred)
    
print("Mean Squared Error:", mse)

winsound.PlaySound(directory+"//#dev//situation_log_updated.wav", winsound.SND_FILENAME)

#joblib.dump(rsf, directory+"//random_forest_survivor.joblib")

plt.scatter(y_actual, y_pred)
plt.xlabel('Actual Values')
plt.ylabel('Predicted Values')
plt.title('Actual vs Predicted')
plt.plot([y_actual.min(), y_actual.max()], [y_actual.min(), y_actual.max()], 'k--')
plt.show()