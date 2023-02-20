#bash script to submit for all 480 combo run

source ~/.bashrc

#module add perl/5.18.2
module add bowtie/1.2.0
module add python/2.7.12

export PERL5LIB=$PERL5LIB:/nas/longleaf/home/mikelaff/tools/mirdeep2_0_0_8/lib
export PERL5LIB=$PERL5LIB:/nas/longleaf/home/mikelaff/perl5/lib/perl5

jid2=$(ssub --notify=ON --mem=64g --wrap=\"miRDeep2.pl combo_reads.fa /proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/data/genome/hg38_UCSC/bowtie_index/hg38.fa combo_v_genome.arf /proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/data/mirbase22/hsa_mature.fa /proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/data/mirbase22/related_hsa_mature.fa /proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/data/mirbase22/hsa_hairpin.fa -t hsa\")
