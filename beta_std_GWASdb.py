import os
import sys
import math
import dateutil
#import pandas as pd
import matplotlib.pyplot as plt
from gwas_db.gwassql import GwasMySqlDb
db = GwasMySqlDb()

#fig = plt.figure()
#fig.set_size_inches(60, 60)
#plot_num = 0
#un_sd_betas = []
#sd_beta_vars = []
#sd_beta_ses = []
#est_vars = []
def extract_single_study(cur, gwas_name):
    cur.execute(
        'SELECT beta_for_alt_allele, betase,af,n_case,n_control,n_total'
        ' FROM gwas_effects JOIN rsid_table AS r USING (rsid)'
        " WHERE Study_ID='{}'".format(gwas_name)
    )
    return cur.fetchall()
fig = plt.figure()
fig.set_size_inches(150, 150)
for i in range(2, 200):
    study = extract_single_study(db.cur, 'GWAS{}'.format(i))
    if study[0][3] == "NULL" and study[0][4] == "NULL":
        continue
    else:
        un_sd_betas = []
        sd_betas = []
        for j in range(0, len(study)):
            if study[j][5] and study[j][5] != "NULL" and study[j][0]:
                un_sd_beta = float(study[j][0])
                se = float(study[j][1])
                freq = float(study[j][2])
                samples = int(study[j][5])
                if se > 0 and samples > 0 and freq > 0:
                    if freq < 1:
                        sd_beta = un_sd_beta / math.sqrt(se * se * 2 * samples * freq * (1 - freq))
                        un_sd_betas.append(un_sd_beta)
                        sd_betas.append(sd_beta)
    if len(un_sd_betas) > 1 and len(un_sd_betas) > 1:
        plt.subplot(15, 15, i)
        plt.hist(un_sd_betas, color="k", alpha=0.5, bins=20)
        plt.hist(sd_betas, color="r", alpha=0.5, bins=20)
fig.savefig('GWASdb.png')






