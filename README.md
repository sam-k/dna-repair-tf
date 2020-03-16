# dna-repair-tf
Genomic analysis of decreased DNA repair at TF-binding sites

Description coming soon!

## Data Sources

### Somatic Mutations: ICGC Data Portal

Somatic mutation data is from [ICGC DCC Data Release 28](https://dcc.icgc.org/releases/release_28), published by the [International Cancer Genome Consortium Data Portal](https://dcc.icgc.org/).

I used the following cancer type datasets:

* BLCA: BLCA-CN (Bladder Cancer - CN), BLCA-US (Bladder Urothelial Cancer - TCGA, US)
* BRCA: BRCA-EU (Breast ER+ and HER2- Cancer - EU/UK), BRCA-FR (Breast Cancer - FR), BRCA-KR (Breast Cancer - Very young women - KR), BRCA-UK (Breast Triple Negative/Lobular Cancer - UK), BRCA-US (Breast Cancer - TCGA, US)
* COAD: COAD-US (Colon Adenocarcinoma - TCGA, US)
* COCA: COCA-CN (Colorectal Cancer - CN)
* HNSC: HNSC-US (Head and Neck Squamous Cell Carcinoma - TCGA, US)
* LUAD: LUAD-US (Lung Adenocarcinoma - TCGA, US)
* LUSC: LUSC-CN (Lung Cancer - CN), LUSC-KR (Lung Cancer - KR), LUSC-US (Lung Squamous Cell Carcinoma - TCGA, US)
* MELA: MELA-AU (Skin Cancer - AU)
* READ: READ-US (Rectum Adenocarcinoma - TCGA, US)
* SKCA: SKCA-BR (Skin Adenocarcinoma - BR)
* SKCM: SKCM-US (Skin Cutaneous melanoma - TCGA, US)

### Active TFBSs: BBGLab

Active transcription factor-binding site data is that used by the [Barcelona Biomedical Genomics Lab](http://bg.upf.edu/group/projects/tfbs/) in their [Sabarinathan et al., 2016](https://www.doi.org/10.1038/nature17661) study.

I paired the following active TFBS (DHS) datasets with the following somatic mutation datasets:

* Proximal BRCA (Breast): BRCA
* Proximal BLCA (Bladder): BLCA
* Proximal CRC (Colorectal): COAD, COCA, READ
* Proximal HNSC (Head and neck squamous cell carcinoma): HNSC
* Proximal LUAD/LUSC (Lung adenocarcinoma, lung squamous cell carcinoma): LUAD, LUSC
* Proximal SKCM (Melanoma): MELA, SKCA, SKCM
* Distal SKCM (Melanoma): MELA, SKCA, SKCM
