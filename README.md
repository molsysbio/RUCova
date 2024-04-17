## RUCova: Removal of Unwanted Covariance in mass cytometry data
Here we present the R package RUCova, a novel method designed to address confounding factors such as heterogeneous cell size and staining efficiency in mass cytometry data. RUCova  removes unwanted covariance using multivariate linear regression based on Surrogates of Unwanted Covariance (SUCs), and Principal Component Analysis (PCA). 

### 1. Install RUCova

Simply run the following in R:

```
remotes::install_github("molsysbio/RUCova")
library(RUCova)
```

### 2. Define the data frame, SUCs, and markers

```data``` should be a data frame containing the single-cell marker signals [rows = cells, columns = markers, features, metadata].


### 3. Apply RUCova

### 4. Evaluate the benefit of RUCova


