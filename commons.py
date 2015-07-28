import sys
import gzip

tot_KG = 0
tot_UKBIOBANK = 0
tot_ALSPAC = 0
tot_common = 0
tot_gwascat = 0
output_file = open(sys.argv[1], "a")

for i in range(1, 23):
    print i
    UKBIOBANK = {}
    KG = {}
    ALSPAC = {}
    GWASCAT = {}
    TWINSUK = {}
    for line in gzip.open("/data/library/1000Genomes/IMPUTE2/1000GP_Phase3/1000GP_Phase3_chr" + str(i) + ".legend.gz",
                          "rb"):
        line = ' '.join(line.split())
        line = line.split(" ")
        if line[0].split(":")[0].startswith('rs'):
            KG[line[0].split(":")[0]] = line[1] + ":" + line[2] + ":" + line[3]

    for line in gzip.open(
            '/data/scratch/project/pharma/ukbiobank/imputed_association_results/asthma/asthma_chr{0}'
            '.assoc.logistic.QCed.gz'.format(str(i)), "rb"):
        line = ' '.join(line.split())
        line = line.split(" ")
        UKBIOBANK[line[1]] = line[4] + ":" + line[2] + ":" + line[3]  #rsid

    for line in open("/data/scratch/project/uk10k/alspac_VCFs_frq_0.01.txt", "rb"):
        line = ' '.join(line.split())
        line = line.split(" ")
        if line[0] == str(i):
            ALSPAC[line[1]] = line[1]

    for line in open("/data/library/GWASCatalog/gwas_catalog_v1.0.1-downloaded_2015-06-22.tsv", "rb"):
        line2 = line.strip()
        line2 = line2.split("\t")
        if line2[21] != "" and line2[11] == str(i):
            GWASCAT[line2[21]] = line
    common = 0
    for pos in GWASCAT:
        print pos
        if pos in KG and pos in UKBIOBANK and KG[pos].split(":")[0] in ALSPAC:
            common += 1
        else:
            print >>output_file, GWASCAT[pos].rstrip('\n')
    print len(GWASCAT), common