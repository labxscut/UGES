# The data for UGES

The original data are available in TCGA (https://gdc.cancer.gov/about-data/publications/pancanatlas):

· ① DNA Methylation (450K Only) - jhu-usc.edu_PANCAN_HumanMethylation450.betaValue_whitelisted.tsv,

and two public datasets of cBioportal platform (http://www.cbioportal.org/):

· ② Breast Cancer (METABRIC, Nature 2012 & Nat Commun 2016)(https://www.cbioportal.org/study/summary?id=brca_metabric)

· ③ Breast Invasive Carcinoma (TCGA, PanCancer Atlas)(https://www.cbioportal.org/study/summary?id=brca_tcga_pan_can_atlas_2018)

The datasets downloaded from cbioportal and the pre-processing data are restored in <https://github.com/labxscut/UGES> or in the latest release named UGES at <https://github.com/labxscut/UGES>. The readers can download them directly.

## Brief description for the data and their potential use

The dataset ② and ③ are two large-scale multi-omics datasets of breast cancer patients, containing the genetic, epigenetic, expression and clinical data, etc. The dataset ① actually is the DNA methylation data of samples included in dataset ③, as the supplementary data. 

In UGES, we used these large-scale multi-omics data as the original data. <Combined_data.csv> and <Standard_combined_data.csv> are the pro-precessing data. <survival.csv> and <UGES_survival.csv> are data used in survival analysis. They provideed us the sufficient condition to do the analysis.

## Data types

In UGES, the data are:

· The alterations of multi-omics data, including mutation, copy number abberation, methylation alterations, etc.

· The clinical data of samples, including intrinsic subtype, age, sex, patient ID, etc.

## Estimate of dataset size

After pre-processing, we have a few GB size of data for 2065 breast cancer patients.
