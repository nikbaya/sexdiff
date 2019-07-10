#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 13:43:33 2019

Identify SNPs that are significantly different for males vs. females for all
phenotypes in rg_tb.

@author: nbaya
"""

import hail as hl
import datetime as dt
import subprocess
import argparse

parser = argparse.ArgumentParser(add_help=False)
parser.add_argument('--parsplit', type=int, required=True, help="number of batches to split phsource phenotypes into")
parser.add_argument('--paridx', type=int, required=True, help="batch id number")

args = parser.parse_args()

parsplit = args.parsplit
paridx = args.paridx

path = 'gs://ukb-mega-gwas-results/round2/additive-tsvs/'
wd = 'gs://nbaya/sexdiff/'

rg_tb = hl.import_table('gs://nbaya/ukb31063.rg_sex.gwas_v2.h2_z4.18jun2019.tsv',impute=True)

phens = rg_tb.phenotype.collect()

start_idx = int((paridx-1)*len(phens)/args.parsplit)
stop_idx = int((paridx)*len(phens)/parsplit)
idx = range(start_idx,stop_idx,1) #chunks all phenotypes for phsource into parsplit number of chunks, then runs on the paridx-th chunk

for phen_i in idx:
    phen = phens[phen_i]
    output_path = wd+f'{phen}.diffgwasloci.tsv.bgz'
    try:
        subprocess.check_output([f'gsutil','ls',output_path]) != None
        print(f'\n#############\n{phen} already completed!\n#############\n')
    except:
        print(f'\n#############\nStarting phenotype {phen} ({idx.index(phen_i)+1} of {len(phens)} for paridx {paridx})\n#############\n')
        start = dt.datetime.now()
        f = hl.import_table(path+phen+'.gwas.imputed_v3.female.tsv.bgz',force_bgz=True,impute=True,key='variant')
        m = hl.import_table(path+phen+'.gwas.imputed_v3.male.tsv.bgz',force_bgz=True,impute=True,key='variant')
        both = f.join(m)
        both1 = both.filter(~(both.low_confidence_variant | both.low_confidence_variant_1)) #remove low confidence variants
        both2 = both1.annotate(diff = (both1.beta-both1.beta_1),
                               diff_se = hl.sqrt(both1.se**2 + both1.se_1**2))
        both3 = both2.annotate(diff_pval = 2*hl.pnorm(-hl.abs(both2.diff/both2.diff_se)))
        both3.select('diff','diff_se','diff_pval').export(output_path)
        print(f'\n#############\nTime for phenotype {phen}: {round((dt.datetime.now()-start).seconds/60, 2)} min\n#############')