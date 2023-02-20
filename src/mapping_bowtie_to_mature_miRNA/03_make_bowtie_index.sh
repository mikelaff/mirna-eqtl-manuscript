#!/bin/bash


fasta_reference="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/data/mirge2.0/miRge.Libs/human/fasta.Libs/human_mirna_SNP_pseudo_miRBase.fa"
output_base="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/data/bowtie_index/miRge_human_mirna_SNP_pseudo_miRBase/miRge"

source ~/.bashrc
module add bowtie/1.2.3

ssub --mem=16g -p steinlab --wrap=\"bowtie-build -f ${fasta_reference} ${output_base}\"

