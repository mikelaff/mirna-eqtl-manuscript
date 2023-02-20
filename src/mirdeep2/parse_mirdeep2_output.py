import os
import csv

mirdeep_output_file = "/Users/mikelaff/Projects/smallRNA/mirdeep2/combo_mirdeep2/mirdeep_runs/run_28_06_2017_t_17_54_34/output.mrd"
working_dir = "/Users/mikelaff/Projects/smallRNA/mirdeep2/combo_mirdeep2/mirdeep_runs/run_28_06_2017_t_17_54_34/"

os.chdir(working_dir)
print("Current working directory: {}".format(os.getcwd()))

#open file for reading
with open(mirdeep_output_file, 'r') as f:
    output_file = f.read().splitlines()

ids = []
score_tot = []
score_star = []
score_reads = []
score_mfe = []
score_rand = []
score_cons = []
seed = []
reads_tot = []
reads_mat = []
reads_loop = []
reads_star = []

for line in range(0,len(output_file)):
    # if line starts with >
    if output_file[line].startswith(">"):
        # get id
        ids.append(output_file[line][1:])
        # get total score
        score_tot.append(output_file[line+1].split()[-1])
        # get star score
        score_star.append(output_file[line + 2].split()[-1])
        # get reads score
        score_reads.append(output_file[line + 3].split()[-1])
        # get mfe score
        score_mfe.append(output_file[line + 4].split()[-1])
        # get rand score
        score_rand.append(output_file[line + 5].split()[-1])
        # get cons score
        score_cons.append(output_file[line + 6].split()[-1])
        # get seed
        if output_file[line + 7].startswith("miRNA"):
            bump = 0
            seed.append(output_file[line + 7].split()[-1])
        else:
            bump = 1
            seed.append("NA")
        # get reads tot
        reads_tot.append(output_file[line + 8 - bump].split()[-1])
        # get reads mat
        reads_mat.append(output_file[line + 9 - bump].split()[-1])
        # get reads loop
        reads_loop.append(output_file[line + 10 - bump].split()[-1])
        # get reads star
        reads_star.append(output_file[line + 11 - bump].split()[-1])

print("Length of ids: {}".format(len(ids)))
print("Length of score_tot: {}".format(len(score_tot)))
print("Length of score_star: {}".format(len(score_star)))
print("Length of score_reads: {}".format(len(score_reads)))
print("Length of score_mfe: {}".format(len(score_mfe)))
print("Length of score_rand: {}".format(len(score_rand)))
print("Length of score_cons: {}".format(len(score_cons)))
print("Length of seed: {}".format(len(seed)))
print("Length of reads_tot: {}".format(len(reads_tot)))
print("Length of reads_mat: {}".format(len(reads_mat)))
print("Length of reads_loop: {}".format(len(reads_loop)))
print("Length of reads_star: {}".format(len(reads_star)))

print(ids[29:32])
print(score_tot[29:32])
print(score_star[29:32])
print(score_reads[29:32])
print(score_mfe[29:32])
print(score_rand[29:32])
print(score_cons[29:32])
print(seed[29:32])
print(reads_tot[29:32])
print(reads_mat[29:32])
print(reads_loop[29:32])
print(reads_star[29:32])

# write csv
with open("mirdeep_scores_parsed.csv", 'w') as output:
    wr = csv.writer(output, dialect='excel')
    headers = ["id", "score_tot", "score_star", "score_reads", "score_mfe", "score_rand", "score_cons", "seed", "reads_tot", "reads_mat", "reads_loop", "reads_star"]
    wr.writerows([headers])
    rows = zip(ids, score_tot, score_star, score_reads, score_mfe, score_rand, score_cons, seed, reads_tot, reads_mat, reads_loop, reads_star)
    for row in rows:
        wr.writerows([row])