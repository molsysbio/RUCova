#' HNSCC data set
#'
#' Is a tibble containing mass cytometry data of single-cell marker signals (rows = cells, columns = markers and metadata) in linear scale. 
#' This data set should be clean, meaning you excluded beads, debris, doublets, dead cells, and single-cells are demultiplexed (important if you want to adapt the linear fits to the samples).
#' In this example we offer a mass cytometry data set consisting of 8 Head-and-Neck Squamous Cell Carcinoma (HNSCC) lines in irradiated (10 Gy) and control (0 Gy) conditions (Figure 2 and Figure 3 in the manuscript).
#'
#' @docType data
#' @format A data frame with 108649 rows and 59 variables:
#' #' \describe{
#'   \item{CXCL1}{Signal for the CXCL1 marker, a cytokine associated with inflammatory responses.}
#'   \item{Ce140Di}{Signal for the Cerium 140 marker, typically used as a control or calibration marker.}
#'   \item{Cyclin_B1}{Marker for the G2/M cell cycle phase regulator Cyclin B1.}
#'   \item{Cyclin_D1}{Marker for the G1 cell cycle phase regulator Cyclin D1.}
#'   \item{DNA_191Ir, DNA_193Ir}{Signals for DNA-intercalating markers labeled with Iridium isotopes, used to assess nuclear content.}
#'   \item{Dead_cells_194Pt, Dead_cells_195Pt, Dead_cells_196Pt, Dead_cells_198Pt}{Signals for dead cell markers labeled with Platinum isotopes, used to identify and exclude dead cells.}
#'   \item{Event_length}{Measure of event duration in the mass cytometer, used for quality control.}
#'   \item{GDF15}{Signal for GDF15, a marker associated with stress and inflammation.}
#'   \item{IdU}{Signal for IdU (Iododeoxyuridine), used to assess DNA synthesis.}
#'   \item{IkBa}{Signal for IkBa, an inhibitor of the NF-kB signaling pathway.}
#'   \item{Ki.67}{Signal for Ki-67, a marker of cell proliferation.}
#'   \item{Lamin_B1}{Signal for Lamin B1, a nuclear lamina protein.}
#'   \item{NICD}{Signal for the Notch Intracellular Domain, a marker of Notch signaling activity.}
#'   \item{Pd102Di, Pd104Di, Pd105Di, Pd106Di, Pd108Di, Pd110Di}{Signals for Palladium isotopes, used as barcodes for cell multiplexing.}
#'   \item{Pt194Di_norm}{Normalized signal for Platinum 194 isotope.}
#'   \item{Time}{Timestamp for the acquisition of each event.}
#'   \item{YAP}{Signal for Yes-associated protein (YAP), a transcriptional co-activator in the Hippo pathway.}
#'   \item{bc}{Barcode channel used for multiplexing.}
#'   \item{beadsOut}{Indicator for exclusion due to bead contamination.}
#'   \item{cCasp3}{Signal for cleaved Caspase-3, a marker of apoptosis.}
#'   \item{cPARP}{Signal for cleaved PARP, another marker of apoptosis.}
#'   \item{dose}{Treatment condition: either \code{0Gy} (control) or \code{10Gy} (irradiated).}
#'   \item{irradiated}{Indicator for whether the sample was irradiated.}
#'   \item{line}{HNSCC cell lines included in the data: "Cal27", "Cal33", "UPCISCC099", "UPCISCC131", "UTSCC16A", "UDSCC2", "UPCISCC154", and "VUSCC147".}
#'   \item{lowPt}{Indicator of low Platinum signal.}
#'   \item{p.p38, p.p53, p4e.BP1}{Signals for phosphorylated p38, p53, and 4E-BP1, markers of stress response and translation regulation.}
#'   \item{pAkt, pAkt_T308}{Signals for phosphorylated Akt, a marker of PI3K/Akt pathway activity.}
#'   \item{pCDC25c, pChk2}{Signals for phosphorylated CDC25c and Chk2, markers of DNA damage response.}
#'   \item{pERK1.2}{Signal for phosphorylated ERK1/2, a marker of MAPK pathway activity.}
#'   \item{pH2A.X, pH3}{Signals for phosphorylated H2A.X and H3, markers of DNA damage and mitosis, respectively.}
#'   \item{pMEK1.2, pNFkB}{Signals for phosphorylated MEK1/2 and NF-kB, markers of MAPK and inflammatory signaling.}
#'   \item{pRb}{Signal for phosphorylated Rb, a marker of cell cycle regulation.}
#'   \item{pS6}{Signal for phosphorylated S6, a marker of protein synthesis.}
#'   \item{pSmad1.8, pSmad2.3}{Signals for phosphorylated Smad1/8 and Smad2/3, markers of TGF-beta signaling.}
#'   \item{pStat1, pStat3}{Signals for phosphorylated Stat1 and Stat3, markers of JAK/STAT pathway activity.}
#'   \item{pan_Akt}{Total Akt protein signal.}
#'   \item{replicate}{Identifier for technical replicates.}
#'   \item{singlets}{Indicator for singlet events, excluding doublets.}
#'   \item{total_ERK}{Total ERK protein signal.}
#' }
"HNSCC_data"