# scRNAseq-periodontitis-analysis
CSN Lab, SNU

Minimal template scripts outlining the scRNA-seq analysis workflow used in the associated periodontitis study.  
These scripts are non-executable examples and do not contain raw data, full code, or patient-related information.

# Required R environment for scRNAseq-periodontitis-analysis
R >= 4.3.3

Seurat==5.2.1
Harmony==1.2.3
SingleR==2.4.1
celldex==1.12.0
Monocle3==1.3.7
ggplot2==3.4.4

# Common dependencies
Matrix
dplyr
tidyr
patchwork
RANN
uwot
BiocManager


## Repository Structure
```plaintext
scRNAseq-periodontitis-analysis/
├── 01_preprocessing.R
├── 02_Integration.R
├── 03_Iterative_Clustering.R
└── 04_Trajectory.R
