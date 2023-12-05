# Predicting Patient Survival

* 119 Patients without TP53 mutations were selected from the Ovarian Serous Cystadenocarcinoma Dataset (PanCanAtlas, 2018)
* The data used included mRNA expression, data on the presence of mutations, and the levels of p53 phosphorylation obtained via the ODE model.
* The Cox proportional-hazards model was used to select features important for patient survival.
* The Random Forest Survivor model with Bayes search optimizer was trained on synthetic data obtained from 102/119 patients from the original data set, as 119 patents on its own were not enough to train the predictor. The model was tested on the remaining 17 patients. Additionally, all the patients were stratified by their overall survival to maintain the balance of classes.
* Mean Squared Error (MSE) metric was chosen for estimating the prediction of the trained model. The best-trained model has been able to reduce MSE error to the value of 355.

![image](https://github.com/Aetherum17/AI-for-Personalised-Medicine/assets/46795020/23937bbe-af4e-431c-b3ce-be27386d48c4)

