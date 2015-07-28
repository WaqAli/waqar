import os
import sys
import math
import pandas as pd
import matplotlib.pyplot as plt


def read_sample_file(sample_file):
    return pd.read_csv(sample_file, sep="\t")

freqs = {}
for line in open('/data/scratch/project/uk10k/twinsuk_VCFs_frq_0.01.txt', 'r'):
    line = line.strip().split()
    if line[0] == '1':
        freqs[line[1]] = line[4]

samples = read_sample_file(sys.argv[1])


fig = plt.figure()
fig.set_size_inches(60, 60)
plot_num = 0
for pheno in samples.columns.values[15:65]:
    plot_num += 1
    pheno_var = samples[pheno].var()
    print pheno, pheno_var
    un_sd_betas = []
    sd_beta_vars = []
    sd_beta_ses = []
    est_vars = []
    if os.path.isfile('/data/scratch/project/uk10k/association_result_twinsuk/1.TWINSUK.beagle.anno.csq.shapeit'
        '.20131101_0.05.{}.assoc.linear'.format(pheno)):
        assoc_file = open('/data/scratch/project/uk10k/association_result_twinsuk/1.TWINSUK.beagle.anno.csq.shapeit'
        '.20131101_0.05.{}.assoc.linear'.format(pheno),"r")
        for line in assoc_file:
            line = line.strip().split()
            if line[2] in freqs and line[4] == "ADD":
                un_sd_beta = float(line[6])
                se = un_sd_beta / float(line[7])
                freq = float(freqs[line[2]])
                num_samples = int(line[5])
                sd_beta_var = un_sd_beta / math.sqrt(pheno_var)
                est_var = se * se * (2 * num_samples * freq * (1 - freq))
                est_vars.append(est_var)
                sd_beta_se = un_sd_beta / math.sqrt(est_var)
                un_sd_betas.append(un_sd_beta)
                sd_beta_vars.append(sd_beta_var)
                sd_beta_ses.append(sd_beta_se)
        plt.subplot(8, 8, plot_num)
 #   plt.axvline(pheno_var, color='r')
        plt.hist(un_sd_betas, color="k", alpha=0.5, bins=20)
        plt.hist(sd_beta_vars, color="r", alpha=0.5, bins=20)
        plt.hist(sd_beta_ses, color="b", alpha=0.5, bins=20)
        plt.xlim((-0.3, 0.3))
        assoc_file.close()

  #  plt.hist(sd_beta_vars, color="r", alpha=0.5)
 #   plt.hist(sd_beta_ses, color="b", alpha=0.5)
fig.savefig('foo.png')






