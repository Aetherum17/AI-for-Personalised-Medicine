import os
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
from sklearn.ensemble import RandomForestRegressor
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error
import numpy as np
import winsound

from sdv.single_table import CTGANSynthesizer
from sdv.metadata import SingleTableMetadata

# Setting up the Paths ############################################################################
directory = os.getcwd()

generate_synthetic_data = False
###################################################################################################

# Load the dataset
data = pd.read_csv(directory+"\\output_file.csv", sep = ";", decimal=',')
data = data[data['SURV_STATUS'] == 1]
data.drop(["SAMPLE_ID", "SURV_STATUS"], axis=1, inplace=True)

# Generate Synthetic Data #########################################################################

# Split the data for generating synthetic data and testing
data_train, data_test = train_test_split(data, test_size=0.10, random_state=42)

if generate_synthetic_data == True:
    # Create metadata from the data for synthetic generation
    metadata = SingleTableMetadata()
    metadata.detect_from_dataframe(data=data_train)
    print("Obtained Metadata")
    
    # Define Synthesizer
    synthesizer = CTGANSynthesizer(metadata, enforce_rounding=False, enforce_min_max_values = False, verbose = False, epochs = 100)
    # Fit Synthesizer
    synthesizer.fit(data_train)
    # Obtain Synthetic data
    synthetic_data = synthesizer.sample(num_rows=5000)
    
    print(synthetic_data.size)
    synthetic_data.to_csv('synthetic_data.csv', index=False)
else:
    synthetic_data = pd.read_csv(directory+"\\synthetic_data_73.csv", sep = ",", decimal='.')
    
# split synthetic and testing data into X and y 
synthetic_data_X = synthetic_data.drop("OS_MONTHS", axis=1).copy()
synthetic_data_y = synthetic_data['OS_MONTHS']
data_test_X = data_test.drop("OS_MONTHS", axis=1).copy()
data_test_y = data_test['OS_MONTHS']

# Define Random Forest Regressor
rf = RandomForestRegressor()
rf.fit(synthetic_data_X, synthetic_data_y)

# Important Features #############################################################################

# Select Important features
importances = rf.feature_importances_
sorted_indices = importances.argsort()[::-1]
feature_names = synthetic_data_X.columns
print(importances)

# Plot the Important features
plt.figure(figsize=(10, 6))
plt.bar(range(len(importances)), importances[sorted_indices], align='center')
plt.xticks(range(len(importances)), feature_names[sorted_indices], rotation=90)
plt.xlabel('Features')
plt.ylabel('Importance')
plt.title('Feature Importance')
plt.tight_layout()
plt.show()

# Set the threshold to select only important features
threshold = 0.01
important_indices = sorted_indices[importances[sorted_indices] > threshold]
important_features = feature_names[important_indices]
synthetic_data_filtered = synthetic_data_X[important_features]

data_test_X = data_test_X[important_features]
y_test = data_test_y

print(synthetic_data_filtered.columns)

# Grid search to define the best parameters
def perform_grid_search():
    param_grid = {
        #'max_depth': [1, 2, 3, 5, 7, 10],
        #'min_samples_split': [1, 2, 3, 5, 7, 10],
        #'min_samples_leaf': [1, 2, 3, 5, 7, 10]
        #'min_samples_leaf': [10, 12, 15, 17, 20, 21]
    }
    grid_search = GridSearchCV(estimator=rf, param_grid=param_grid, scoring='explained_variance', cv=5)
    grid_search.fit(synthetic_data_X, synthetic_data_y)
    
    best_params = grid_search.best_params_
    best_score = grid_search.best_score_
    
    print("Best Parameters:", best_params)
    print(best_score)
    winsound.PlaySound(directory+"//#dev//situation_log_updated.wav", winsound.SND_FILENAME)

def build_model():
    rf = RandomForestRegressor(max_depth = 7, min_samples_split = 3, min_samples_leaf = 10, random_state = 42)
    
    rf.fit(synthetic_data_filtered, synthetic_data_y)
    print("Model Fit")
    use_model_for_prediction(model = rf)

def use_model_for_prediction(model):
    y_pred = model.predict(data_test_X)
    #print("Predicted Values")
    
    mse = mean_squared_error(y_test, y_pred)
    
    print("Mean Squared Error:", mse)
    
    plt.scatter(y_test, y_pred)
    plt.xlabel('Actual Values')
    plt.ylabel('Predicted Values')
    plt.title('Actual vs Predicted')
    plt.plot([y_test.min(), y_test.max()], [y_test.min(), y_test.max()], 'k--')
    plt.show()
    
    winsound.PlaySound(directory+"//#dev//situation_log_updated.wav", winsound.SND_FILENAME)
    
#perform_grid_search()
build_model()