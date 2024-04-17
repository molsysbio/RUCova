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

#### mean BC:

The RUCova function called ```RUCova::calc_mean_BC``` works in two steps. First, it applies the ```asinh()``` function and adjusts the transformed distributions of the barcoding isotopes by matching a specific percentile ```q```. Then, it looks at the signals of each isotope for each cell and picks the top ```n_bc``` signals used for barcoding. The function returns a vector with the mean BC signals in linear scale, as it applyies the inverse tranformation ```sinh()```. In the following example the Cell-ID 20-Plex Pd Barcoding Kit was used, where 3 out of 6 isotopes are mixed together to form one barcode:
```
 data |> 
  mutate(mean_BC = RUCova::calc_mean_BC(Pd102Di, Pd104Di, Pd105Di, Pd106Di, Pd108Di, Pd110Di,
                                        n_bc = 3, q = 0.95))
```

#### mean DNA:

The RUCova function called ```RUCova::calc_mean_DNA``` applies the ```asinh()``` function and adjusts the transformed distributions of the iridium isotopes (Ir191 and Ir193) by matching a specific percentile ```q``` to then take the mean value of these two signals per cell. The function returns a vector with the mean DNA signals in linear scale, as it applyies the inverse tranformation ```sinh()```.
```
 data |> 
  mutate(mean_DNA = RUCova::calc_mean_DNA(DNA_191Ir, DNA_193Ir, q = 0.95)
```

Then we define the vector for the SUCs: ```surrogates = c("total_ERK", "pan_Akt", "mean_DNA", "mean_BC")```.

SUCs can also be linearly transformed onto a new coordinate system by applying Principal Component Analysis. By doing this, users can decide the extent of unwanted correlation to be removed by using eg.: only PC1 as a predictive variable in the univariate model, or PC1 to PC2, or all PCs which is equivalent as taking all surrogates as predictive variables.

#### PCA on SUCs:

```
pca_sucs <- data |> 
  select(surrogates) |> 
  mutate_all(asinh) |> 
  mutate_all(scale) |> 
  prcomp()
```

Calculate and plot the variance explained by each PC:

```
tibble(perc = as.numeric(pca_sucs$sdev^2/sum(pca_sucs$sdev^2))*100, 
       PC = 1:length(pca_sucs$sdev)) |>  
  ggplot(aes(x = PC, y = perc, fill = as.character(PC), label = round(perc,1))) + 
  geom_col() +
  geom_label()
```

Check the loadings of each PC:

```
as.data.frame(pca_sucs$rotation) |> 
  rownames_to_column("SUC") |> 
  pivot_longer(names_to = "PC", values_to = "loadings", -SUC) |> 
  ggplot(aes(x = loadings, y = SUC)) + 
  geom_col() + 
  facet_wrap(~PC, nrow = 1)
```

Add the PCs to the main data frame:

```
data <- data |> 
        cbind(as.data.frame(pca_sucs$x)) 
```

Finally we define a vector for the markers for which we want to regress-out the unwated covariance: ```markers = c("pH3","IdU","Cyclin_D1","Cyclin_B1", "Ki.67","pRb","pH2A.X","p.p53","p.p38","pChk2","pCDC25c","cCasp3","cPARP","pAkt","pAkt_T308","pMEK1.2","pERK1.2","pS6","p4e.BP1","pSmad1.8","pSmad2.3","pNF.κB","IκBα", "CXCL1","Lamin_B1", "pStat1","pStat3", "YAP","NICD")```
The character variables in the vectors ```surrogates```and ```markers``` must be findable as column names in the mass cytometry data set ```data```.

### 3. Apply RUCova

Let's imagine you want to be conservative and only remove correlations between markers and PC1 (of SUCs). Then, ```SUCs= "PC1"```,  ```apply_asinh_SUCs = FALSE``` as asinh transformation is not necessary on PCs (it was applied on SUCs before PCA). We will choose the interaction model (```model = "interaction"```), as the data set contains different cell lines (samples), and for each one we want to allow different slopes and intercepts between marker expression and SUCs (or PCs (```col_name_sample = "line"```). Differences in marker expression between cell lines can be artificially influenced by e.g.: different cell volumes leading to an unwanted covariance. To remove these artefactual differences in marker expression we center the surrogates across samples for the linear fit (```center_SUCs = "across_samples"```). In case is desirable to keep the remaining offset between the cell lines, we set ``` keep_offset = TRUE```. 

```
rucova_pc1 <-  RUCova::rucova(data, 
                              markers,
                              SUCs= "PC1",
                              apply_asinh_SUCs = FALSE,
                              model = "interaction",
                              col_name_sample = "line",
                              center_SUCs = "across_samples",
                              keep_offset = TRUE)
data_reg_pc1 <- rucova_pc1$data_reg
```

To remove the covariance driven by all PCs (in this example, 4):

```
rucova_all <-  RUCova::rucova(data, 
                              markers,
                              SUCs= c("PC1","PC2","PC3","PC4"),
                              apply_asinh_SUCs = FALSE,
                              model = "interaction",
                              col_name_sample = "line",
                              center_SUCs = "across_samples",
                              keep_offset = TRUE)
data_reg_all <- rucova_all$data_reg
```

More information about the models can be found in: ------------ paper?


### 4. Evaluate the benefit of RUCova

We recommend to calculate the Pearson correlation coefficients between markers and SUCs before RUCova and after RUCova, and compare:

```
corr_reg_before <- data |>
  mutate_at(vars(markers,surrogates), asinh) |> 
  select(markers,surrogates,PC1,PC2,PC3,PC4) |> 
  cor(method= "pearson")

corr_reg_all <- data_reg_all |>
  mutate_at(vars(markers,surrogates), asinh) |> 
  select(markers,surrogates,PC1,PC2,PC3,PC4) |> 
  cor(method= "pearson")
```
Here, ```corr_reg_before```and ```corr_reg_all``` should be square matrices containing the pearson correlation coefficient between markers and surrogates. Row and column names of the matrices should be the corresponding markers and surrogates. Correlation coefficients can be calculated across the entire data set or filtered by each sample if desired.

Using the RUCova function ```heatmap_compare_corr()``` we can generate a square heatmap with the correlation coefficients before RUCova in the lower triangle (first argument) and after RUCova in the upper triangle (second argument) for a direct comparison:

```
RUCova::heatmap_compare_corr(lower = corr_reg_before, upper = corr_reg_all)
```
To further evaluate the benefit of RUCova, we recommend to also perform your favorite analysis on the data before RUCova and also after RUCova, and then compare both. This can be: density plots of marker intensity signals, UMAPs, heatmap of fold-changes, Louvain clustering, etc. 



