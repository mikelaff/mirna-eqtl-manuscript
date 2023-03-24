# MicroRNA-eQTLs in the Developing Human Neocortex

Scripts used in the analyses in, "MicroRNA-eQTLs in the developing human neocortex link miR-4707-3p expression to brain size."  

## Citation

Michael J Lafferty, Nil Aygün, Niyanta K Patel, Oleh Krupa, Dan Liang, Justin M Wolter, Daniel H Geschwind, Luis de la Torre-Ubieta, Jason L Stein (2023) **MicroRNA-eQTLs in the developing human neocortex link miR-4707-3p expression to brain size** eLife 12:e79488.  

https://doi.org/10.7554/eLife.79488  

## Data

### Summary Statistics
Summary statistics for all SNP-miRNA association tests are located at [`results/sumstats/`](https://github.com/mikelaff/mirna-eqtl-manuscript/tree/main/results/sumstats).  

### Small RNA-Sequencing
Raw small RNA sequencing data used in this manuscript will be available via dbGaP accession [phs001900.v1.p1](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs003106.v1.p1).  

### Total RNA-Sequencing
Total RNA sequencing data used in this manuscript is available via dbGaP accession [phs001900.v1.p1](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs001900.v1.p1) and was used in the associated manuscript by [Walker et al.](https://doi.org/10.1016/j.cell.2019.09.021).  

### Microarray Expression
Microarray expression data used in this manuscript is available via GEO accession [GSE224188](https://doi.org/10.1016/j.cell.2019.09.021).

## Code

### Local-miRNA-eQTLs

Primary local-miRNA-eQTLs were found using EMMAX using scripts [here](https://github.com/mikelaff/mirna-eqtl-manuscript/tree/main/src/emmax/final_pipeline). Secondary and greater degree eQTLs were found by repeating the pipeline at various significance thresholds using scripts [here](https://github.com/mikelaff/mirna-eqtl-manuscript/tree/main/src/conditionally_independent_eqtls).





