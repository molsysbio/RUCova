## RUCova: Removal of Unwanted Covariance in mass cytometry data
Here we present the R package RUCova, a novel method designed to address confounding factors such as heterogeneous cell size and staining efficiency in mass cytometry data. RUCova  removes unwanted covariance using multivariate linear regression based on Surrogates of Unwanted Covariance (SUCs), and Principal Component Analysis (PCA). 

### 1. Install RUCova

Simply run the following in R:

```
remotes::install_github("molsysbio/RUCova")
library(RUCova)
```

### 2. Define the data frame, SUCs, and markers

```data```: a data frame containing the single-cell marker signals [rows = cells, columns = markers and metadata] in linear scale. Do not apply any matematical transformation on the signals.

As Surrogates of Unwanted Covariance (SUCs) we selected the signals of total ERK, pan Akt, mean BC (mean value of normalised and used barcoding isotopes), and mean DNA (mean value of iridium DNA intercalators).
Defining mean BC:
The RUCova function called ```RUCova::calc_mean_BC``` works in two steps. First, it adjusts the distributions of the barcoding isotopes, which are used for multiplexing, by matching a specific percentile ```q```. Then, it looks at the signals of each isotope for each cell and picks the top ```n_bc``` signals used for barcoding. In the following example the Cell-ID 20-Plex Pd Barcoding Kit was used, where 3 out of 6 isotopes are mixed together to form one barcode:
```
 data |> 
  mutate(mean_BC = RUCova::calc_mean_BC(Pd102Di, Pd104Di, Pd105Di, Pd106Di, Pd108Di, Pd110Di,
                                        n_bc = 3, q = 0.95))
```






### 3. Apply RUCova

### 4. Evaluate the benefit of RUCova


