#  Short 6-Item Screener for Childbirth-Related PTSD
## Overview

## Data Availability
Due to patient medical confidentiality, we cannot share the raw dataset.

## Running the code
`data_analysis.R` - This R script implements the three steps: 
(1) item selection using L1-penalized logistic regression (LASSO),
(2) model training with bias correction using Firth’s penalized likelihood logistic regression on a train set defined by bootstrap sampling, and 
(3) model evaluation using metrics including Area Under the Receiver Operating Characteristic (ROC) curve (AUC), Sensitivity, Specificity, and Precision (Positive Predictive value: PPV).

`Reduced_PCL.ipynb` - calculate performances of a selecte threshold cutoff for the short 6-item PCL-5 questionnaire
