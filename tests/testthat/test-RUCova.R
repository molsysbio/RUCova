test_that("calc_mean_BC works as expected", {
  sce <- RUCova::sce
  bc_channels <- c("Pd102Di", "Pd104Di", "Pd105Di", "Pd106Di", "Pd108Di", "Pd110Di", "Dead_cells_194Pt", "Dead_cells_198Pt")
  result <- RUCova::calc_mean_BC(sce, name_assay = "counts", bc_channels, n_bc = 4, q = 0.95)
  # Check if the output is a SingleCellExperiment object
  expect_s4_class(result, "SingleCellExperiment")
  # Check for the presence of expected assays
  expect_true("counts" %in% assayNames(result))
  # Check for the presence of mean_BC 
  expect_true("mean_BC" %in% rownames(result))
})


test_that("calc_mean_DNA works as expected", {
  sce <- RUCova::sce
  dna_channels <- c("DNA_191Ir", "DNA_193Ir")
  result <- RUCova::calc_mean_DNA(sce, name_assay = "counts", dna_channels, q = 0.95)
  # Check if the output is a SingleCellExperiment object
  expect_s4_class(result, "SingleCellExperiment")
  # Check for the presence of expected assays
  expect_true("counts" %in% assayNames(result))
  # Check for the presence of mean_BC 
  expect_true("mean_DNA" %in% rownames(result))
})


test_that("rucova function works as expected", {
  sce <- RUCova::sce
  m <- c("pH3","IdU","Cyclin_D1","Cyclin_B1", "Ki.67","pRb","pH2A.X","p.p53","p.p38","pChk2",
          "pCDC25c","cCasp3","cPARP","pAkt","pAkt_T308","pMEK1.2","pERK1.2","pS6","p4e.BP1","pSmad1.8",
          "pSmad2.3","pNFkB","IkBa", "CXCL1","Lamin_B1", "pStat1","pStat3", "YAP","NICD")
  result <- RUCova::rucova(sce = sce, 
                           name_assay_before = "counts", 
                           markers = m, 
                           SUCs = c("total_ERK","pan_Akt"),
                           apply_asinh_SUCs = TRUE,  
                           model = "interaction", 
                           center_SUCs = "across_samples",
                           col_name_sample = "line", 
                           name_assay_after = "counts_interaction")
  # Check if the output is a SingleCellExperiment object
  expect_s4_class(result, "SingleCellExperiment")
  # Check for the presence of expected assays
  expect_true("counts_interaction" %in% assayNames(result))
  # Check for the presence of mean_BC 
})


test_that("heatmap_compare_corr function works as expected", {
  sce <- RUCova::sce
  m <- c("pH3","IdU","Cyclin_D1","Cyclin_B1", "Ki.67","pRb","pH2A.X","p.p53","p.p38","pChk2",
         "pCDC25c","cCasp3","cPARP","pAkt","pAkt_T308","pMEK1.2","pERK1.2","pS6","p4e.BP1","pSmad1.8",
         "pSmad2.3","pNFkB","IkBa", "CXCL1","Lamin_B1", "pStat1","pStat3", "YAP","NICD")
  result <- heatmap_compare_corr(sce[,sce$line == "Cal33"], name_assay_before = "counts")
  # Check if the output is a Heatmap
  expect_s4_class(result, "Heatmap")
})
