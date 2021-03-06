# dna-repair-tf
Genomic analysis of decreased DNA repair at TF-binding sites.

The full report can be found here: [Creating and Measuring Mutation Profiles Within TFBSs](https://www.samkim.io/#/post/tfbs-research-2020)

### Timeline

:checkered_flag: This independent research project was **started** on Jan. 20, 2020, under the guidance of [Dr. Raluca Gordân](https://gordanlab.cs.duke.edu/) and [Harshit Sahay](https://genome.duke.edu/directory/cbb-phd-student/harshit-sahay) at Duke University.

:white_check_mark: This project was **finished** on May 3, 2020.

## Data Sources

### Somatic Mutations: ICGC Data Portal

Somatic mutation data is from [ICGC DCC Data Release 28](https://dcc.icgc.org/releases/release_28), published by the [International Cancer Genome Consortium Data Portal](https://dcc.icgc.org/).

I used the following cancer type datasets:

| Cancer Type | Projects |
|-|-|
| BLCA | BLCA-CN (bladder cancer – CN)<br>BLCA-US (bladder urothelial cancer – TCGA, US) |
| BRCA | BRCA-EU (breast ER+ and HER2- cancer – EU/UK)<br>BRCA-FR (breast cancer – FR)<br>BRCA-KR (breast cancer – very young women – KR)<br>BRCA-UK (breast triple negative/lobular cancer – UK)<br>BRCA-US (breast cancer – TCGA, US) |
| COAD | COAD-US (colon adenocarcinoma – TCGA, US) |
| COCA | COCA-CN (colorectal cancer – CN) |
| HNSC | HNSC-US (head and neck squamous cell carcinoma – TCGA, US) |
| LUAD | LUAD-US (lung adenocarcinoma – TCGA, US) |
| LUSC | LUSC-CN (lung cancer – CN)<br>LUSC-KR (lung cancer – KR)<br>LUSC-US (lung squamous cell carcinoma – TCGA, US) |
| MELA | MELA-AU (skin cancer – AU) |
| READ | READ-US (rectum adenocarcinoma – TCGA, US) |
| SKCA | SKCA-BR (skin adenocarcinoma – BR) |
| SKCM | SKCM-US (skin cutaneous melanoma – TCGA, US) |


### Active TFBSs: BBGLab

Active transcription factor-binding site data is that used by the [Barcelona Biomedical Genomics Lab](http://bg.upf.edu/group/projects/tfbs/) in their [Sabarinathan et al., 2016](https://www.doi.org/10.1038/nature17661) study.

I paired the following active TFBS (DHS) datasets with the following somatic mutation datasets:

| TFBSs | Mutations |
|-|-|
| Proximal BRCA (breast) | BRCA |
| Proximal BLCA (bladder) | BLCA |
| Proximal CRC (colorectal) | COAD, COCA, READ |
| Proximal HNSC (head and neck squamous cell carcinoma) | HNSC |
| Proximal LUAD/LUSC (lung adenocarcinoma, lung squamous cell carcinoma) | LUAD, LUSC |
| Proximal SKCM (melanoma) | MELA, SKCA, SKCM |
| Distal SKCM (melanoma) | MELA, SKCA, SKCM |
