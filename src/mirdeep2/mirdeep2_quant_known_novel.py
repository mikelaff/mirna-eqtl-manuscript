#!/usr/bin/env python3
"""

@author: mikelaff

"""

import os
import subprocess

#globals longleaf
mirdeep_output_root = "/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/mirdeep2/mirbase_and_novel_expression/mirdeep2_runs_by_rnaid/"
fastq_file_list = "/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/data/lists/fastq_files_fullpath.txt"

#globals local
#mirdeep_output_root = "/Users/mikelaff/Repositories/mirna-eqtl/mirdeep2/output/"
#fastq_file_list = "/Users/mikelaff/Repositories/mirna-eqtl/lists/fastq_files_fullpath.txt"

#create output root directory
if not os.path.exists(mirdeep_output_root):
    os.makedirs(mirdeep_output_root)
os.chdir(mirdeep_output_root)
print("Current working directory: {}".format(os.getcwd()))

#load fastq file list
with open(fastq_file_list, 'r') as f:
    fastq_gz_files = f.read().splitlines()
#print("\n".join(fastq_gz_files[0:5]))
#list of fastq files for processing
files=[]
for path in fastq_gz_files:
    files.append(path.split("/")[9])
#print("\n".join(files[0:5]))
#rnaID for directorys
rnaid=[]
for id in files:
    rnaid.append(id.split(".")[0])
#print("\n".join(rnaid[0:5]))

#for each fastq file
for i in range(0,len(files)):
#for i in range(1, 5):
    #restart at output root
    os.chdir(mirdeep_output_root)

    #creat working directory for this new file
    #print("Working on {}".format(files[i]))
    #print("Creating directory: {}".format(rnaid[i]))
    if not os.path.exists("./" + rnaid[i]):
        os.makedirs("./" + rnaid[i])
    #change into this new working directory
    os.chdir("./" + rnaid[i])
    #print("new working dir: {}".format(os.getcwd()))

    #open file to write script
    with open("slurm_{}.sh".format(rnaid[i]), "w") as f:
        f.write("#bash script to submit for file: {}\n".format(rnaid[i]))
        f.write("module add perl/5.18.2\n")
        f.write("module add bowtie/1.2.0\n")
        f.write("module add python/2.7.12\n")
        f.write("source ~/tmp_bashrc\n")

        #convert fastq.gz to fastq
        f.write("jid1=$(ssub --notify=OFF --wrap=\\\"cp {} .\\\")\n".format(fastq_gz_files[i]))
        f.write("jid2=$(ssub --notify=OFF -d afterok:$jid1 --wrap=\\\"gunzip {}\\\")\n".format(files[i]))

        #mapper
        f.write("jid3=$(ssub --notify=OFF -d afterok:$jid2 --mem=8g --wrap=\\\"mapper.pl {}.fastq -e -h -j -k TGGAATTCTCGGGTGCCAAGG -l 18 -m -p /proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/data/genome/hg38_UCSC/bowtie_index/hg38 -s {}_collapsed.fa -t {}_collapsed_vs_genome.arf -v\\\")\n".format(rnaid[i], rnaid[i], rnaid[i]))

        #quantifier
        f.write("jid4=$(ssub --notify=OFF -d afterok:$jid3 --mem=8g --wrap=\\\"quantifier.pl -p /proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/mirdeep2/mirbase_and_novel_expression/fasta_known_and_novel_mirnas/hsa_and_novel_hairpin.fa -m /proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/mirdeep2/mirbase_and_novel_expression/fasta_known_and_novel_mirnas/hsa_and_novel_mature.fa -r {}_collapsed.fa\\\")\n".format(rnaid[i]))

        #mirdeep2
        #f.write("jid5=$(ssub --notify=ON -d afterok:$jid4 --mem=8g --wrap=\\\"miRDeep2.pl {}_collapsed.fa /proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/genome/hg38_UCSC/bowtie_index/hg38.fa {}_collapsed_vs_genome.arf /proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/mirbase/hsa_mature.fa /proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/mirbase/related_hsa_mature.fa /proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/mirbase/hsa_hairpin.fa -t hsa\\\")\n".format(rnaid[i], rnaid[i]))

        #clean up files
        f.write("jid5=$(ssub --notify=OFF -d afterok:$jid4 --wrap=\\\"rm {}.fastq {}_collapsed.fa {}_collapsed_vs_genome.arf\\\")\n".format(rnaid[i], rnaid[i], rnaid[i]))

    print("submitting job for: {}".format(rnaid[i]))
    #subprocess.call(["echo", "ssub --wrap=\\\"sh slurm_{}.sh\\\"".format(rnaid[i])])
    subprocess.call(["sh", "slurm_{}.sh".format(rnaid[i])])

