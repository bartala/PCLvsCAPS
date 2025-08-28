#  Short 6-Item Screener for Childbirth-Related PTSD
## Overview
Childbirth-related posttraumatic stress disorder (CB-PTSD) is an under-recognized maternal mental health condition with significant implications for mothers and infants. While the PTSD Checklist for DSM-5 (PCL-5) is widely validated, its 20 items may be burdensome for postpartum screening. We analyzed data from 107 women who experienced traumatic childbirth, with PTSD diagnosis confirmed via the CAPS-5 gold standard. Using nested bootstrap resampling with LASSO variable selection and Firth-penalized logistic regression, we identified a stable six-item subset of the PCL-5 that retained strong predictive accuracy. The reduced six-item “sum” model achieved AUC = 0.915 (95% CI 0.831–0.980), sensitivity = 0.803, and specificity = 0.816, closely matching the full 20-item scale. A cutoff score of 7 optimized classification (sensitivity = 0.963, specificity = 0.825). This brief screener reduces patient burden and enables rapid triage in postpartum settings, providing a practical tool for early detection of women at high risk for CB-PTSD.

## Data Availability
Due to patient medical confidentiality, we cannot share the raw dataset.

## Running the code
`data_analysis.R` - This R script implements the three steps: 
(1) item selection using L1-penalized logistic regression (LASSO),
(2) model training with bias correction using Firth’s penalized likelihood logistic regression on a train set defined by bootstrap sampling, and 
(3) model evaluation using metrics including Area Under the Receiver Operating Characteristic (ROC) curve (AUC), Sensitivity, Specificity, and Precision (Positive Predictive value: PPV).

`Reduced_PCL.ipynb` - calculate performances of a selecte threshold cutoff for the short 6-item PCL-5 questionnaire

## Miscellaneous
Please send any questions you might have about the code and/or the algorithm to alon.bartal@biu.ac.il.

## Citing
If you find this paper useful for your research, please consider citing us:
```
@article{bartala6item,
  title={Short 6-Item Screener for Childbirth-Related PTSD},
  author={},
  journal={},
  volume={},
  number={},
  pages={},
  year={2025},
  publisher={}
}
```

