## RUCova: Removal of Unwanted Covariance in mass cytometry data
Here we present the R package RUCova, a novel method designed to address confounding factors such as heterogeneous cell size and staining efficiency in mass cytometry data. RUCova  removes unwanted covariance using multivariate linear regression based on Surrogates of Unwanted Covariance (SUCs), and Principal Component Analysis (PCA). 

Citation:

RUCova: Removal of Unwanted Covariance in mass cytometry data
Rosario Astaburuaga-García, Thomas Sell, Samet Mutlu, Anja Sieber, Kirsten Lauber, Nils Blüthgen
Bioinformatics 2024; doi: https://doi.org/10.1101/2024.05.24.595717

### 1. Install RUCova

Run the following in R:

```
remotes::install_github("molsysbio/RUCova@devel", force = TRUE)
library(RUCova)
```

### 2. Define the SingleCellExperiment, SUCs, and markers


We offer a mass cytometry data set consisting of 8 Head-and-Neck Squamous Cell Carcinoma (HNSCC) lines in irradiated (10 Gy) and control (0 Gy) conditions (Figure 2 and Figure 3 in the manuscript) in form of a tibble. This data [rows = markers, columns = cells] is in linear scale and should be clean, meaning without calibration beads, debris, doublets, and dead cells, and single-cells are demultiplexed (important if you want to adapt the linear fits to the samples). 

```
data <- RUCova::HNSCC_data
colnames(data) <- make.names(colnames(data)) 
head(data)
rownames(data) <- 1:dim(data)[1]
data <- data |> mutate(cell_id = 1:n())
```

We convert our data into a ```SingleCellExperiment``` to ensure interoperability across Bioconductor packages. We also add an id per cell to be able to do a cell-wise signals value modification. Metadata of the cells (eg.: samples, treatment, etc) are stored in ``colData``.


```
sce <- SingleCellExperiment(
  assays = list(counts = t(data[,4:48])),    # Assay data
  rowData = DataFrame(measured_channels = colnames(data[,4:48])),# Metadata for rows
  colData = DataFrame(cell_id = data$cell_id,
                      line = data$line,
                      dose = data$dose), # Metadata for columns
)

```

As Surrogates of Unwanted Covariance (SUCs) we selected the signals of total ERK, pan Akt, mean BC (mean value of normalised and used barcoding isotopes), and mean DNA (mean value of iridium DNA intercalators).

#### mean BC

The RUCova function called ```RUCova::calc_mean_BC``` works in two steps. First, it applies the ```asinh()``` function and adjusts the transformed distributions of the barcoding isotopes by matching a specific percentile ```q```. Then, it looks at the signals of each isotope for each cell and picks the top ```n_bc``` signals used for barcoding. The function adds an additional marker (row) with the mean BC signals in linear scale (as it applies the inverse transformation ```sinh()``) in the specified assay. In the following example the Cell-ID 20-Plex Pd Barcoding Kit was used, in addition to Platinum 194 and 198, where 4 out of 8 isotopes are mixed to form one barcode:

```
 bc_channels <- c("Pd102Di", "Pd104Di", "Pd105Di", "Pd106Di", "Pd108Di", "Pd110Di", 
                 "Dead_cells_194Pt", "Dead_cells_198Pt")
sce <- RUCova::calc_mean_BC(sce, name_assay = "counts", bc_channels, n_bc = 4, q = 0.95)
```

#### mean DNA

The RUCova function called ```RUCova::calc_mean_DNA``` applies the ```asinh()``` function and adjusts the transformed distributions of the iridium isotopes (Ir191 and Ir193) by matching a specific percentile ```q``` to take then the mean value of these two signals per cell. The function adds an additional marker (row) with the mean DNA signals in linear scale (as it applies the inverse transformation ```sinh()```) in the specified assay. 

```
dna_channels <- c("DNA_191Ir", "DNA_193Ir")
sce <- RUCova::calc_mean_DNA(sce, name_assay = "counts", dna_channels, q = 0.95)
```

Then we define the vector for the SUCs: ```surrogates = c("total_ERK", "pan_Akt", "mean_DNA", "mean_BC")```.

SUCs can also be linearly transformed into a new coordinate system by applying Principal Component Analysis. By doing this, users can decide the extent of unwanted correlation to be removed by using eg.: only PC1 as a predictive variable in the univariate model, or PC1 to PC2, or all PCs equivalent to taking all surrogates as predictive variables.

#### PCA on SUCs

```
pca <- t(assay(sce,"counts")) |> 
  as.tibble() |> 
  select(surrogates) |> 
  mutate_all(asinh) |> 
  mutate_all(scale) |> 
  prcomp()
```

Calculate and plot the variance explained by each PC:

```
tibble(perc = as.numeric(pca$sdev^2/sum(pca$sdev^2))*100, 
       PC = 1:length(pca$sdev)) |>  
  ggplot(aes(x = PC, y = perc, fill = as.character(PC), label = round(perc,1))) + 
  geom_col() +
  geom_label()
```

Check the loadings of each PC:

```
as.data.frame(pca$rotation) |> 
  rownames_to_column("SUC") |> 
  pivot_longer(names_to = "PC", values_to = "loadings", -SUC) |> 
  ggplot(aes(x = loadings, y = SUC)) + 
  geom_col() + 
  facet_wrap(~PC, nrow = 1)
```

Add the PCs to the SingleCellExperiment object:

```
reducedDim(sce, "PCA") <- pca$x
```

Finally, we define a vector for the markers for which we want to regress out the unwanted covariance: ```m = c("pH3","IdU","Cyclin_D1","Cyclin_B1", "Ki.67","pRb","pH2A.X","p.p53","p.p38","pChk2","pCDC25c","cCasp3","cPARP","pAkt","pAkt_T308","pMEK1.2","pERK1.2","pS6","p4e.BP1","pSmad1.8","pSmad2.3","pNFkB","IkBa", "CXCL1","Lamin_B1", "pStat1","pStat3", "YAP","NICD")```
The character variables in the vectors ```surrogates```and ```m``` must be findable as row names in the mass cytometry assay ```counts```.

### 3. Apply RUCova

Let's imagine you want to be conservative and only remove correlations between markers and PC1 (of SUCs). Then, ```SUCs= "PC1"``` and  ```apply_asinh_SUCs = FALSE```, as asinh transformation is not necessary on PCs (it was applied on SUCs before PCA). We will choose the interaction model (```model = "interaction"```), as the data set contains different cell lines (samples), and for each one, we want to allow different slopes and intercepts between marker expression and SUCs (or PCs (```col_name_sample = "line"```). Differences in marker expression between cell lines can be artificially influenced by e.g.: different cell volumes leading to an unwanted covariance. To remove these artefactual differences in marker expression we center the surrogates across samples for the linear fit (```center_SUCs = "across_samples"```). In case is desirable to keep the remaining offset between the cell lines, we set ``` keep_offset = TRUE```. 

```

sce <- RUCova::rucova(sce = sce, 
               name_assay_before = "counts",
               markers = m,
               SUCs = "PC1",  
               name_reduced_dim = "PCA",
               apply_asinh_SUCs = FALSE, 
               model = "interaction",
               center_SUCs = "across_samples",
               keep_offset = TRUE, 
               col_name_sample = "line",     
               name_assay_after = "counts_interaction")
```

To remove the covariance driven by all PCs (in this example, 4):

```
sce <- RUCova::rucova(sce = sce, 
               name_assay_before = "counts",
               markers = m,
               SUCs =c("PC1","PC2","PC3","PC4"),
               name_reduced_dim = "PCA",
               apply_asinh_SUCs = FALSE, 
               model = "interaction",
               center_SUCs = "across_samples",
               keep_offset = TRUE, 
               col_name_sample = "line",     
               name_assay_after = "counts_interaction_all")
```

The regressed counts are then store in an assay named ``name_assay_after`` and the model details are stored in ``metadata(sce,paste0("model_",name_assay_after))``.

### 4. Evaluate the benefit of RUCova

We recommend calculating the Pearson correlation coefficients between markers and SUCs before RUCova and after RUCova, and comparing. For this we offer a heatmap function that plots the pearson correlation coefficients between markers on a double triangular heatmap (lower triangle: before RUCova, upper triangle: after RUCova). If RUCova has not been applied, the output is a symmetric heatmap.

```
heatmap_compare_corr(sce, name_assay_before = "counts", name_assay_after = "counts_interaction_all", name_reduced_dim = "PCA")
```
Here, ```corr_reg_before```and ```corr_reg_all``` should be square matrices containing the Pearson correlation coefficient between markers and surrogates. Row and column names of the matrices should be the corresponding markers and surrogates. Correlation coefficients can be calculated across the entire data set or filtered by each sample if desired.

Using the RUCova function ```heatmap_compare_corr()``` we can generate a square heatmap with the correlation coefficients before RUCova in the lower triangle (first argument) and after RUCova in the upper triangle (second argument) for a direct comparison:

```
RUCova::heatmap_compare_corr(lower = corr_reg_before, upper = corr_reg_all)
```
To further evaluate the benefit of RUCova, we recommend to also perform your favorite analysis on the data before RUCova and also after RUCova, and then compare both. This can be: density plots of marker intensity signals, UMAPs, heatmap of fold-changes, Louvain clustering, etc. 

For more information about RUCova and examples, check out the vignette:

```
vignette("RUCova")
```

