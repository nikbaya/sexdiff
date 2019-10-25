#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 09:58:45 2019

@author: nbaya
"""

import hail as hl
import pandas as pd

wd = 'gs://nbaya/sexdiff/'

threshold = 5e-8 #significance threshold

rg_tb = hl.import_table('gs://nbaya/ukb31063.rg_sex.gwas_v2.h2_z4.18jun2019.tsv',impute=True)
tb_complete = hl.import_table(wd+'ukb31063.diff_gwas_loci.tsv.bgz',impute=True,force_bgz=True)

phens = rg_tb.phenotype.collect()
descs = rg_tb.description.collect()
n_phens = len(phens)

tb_ls = [None]*n_phens
tb1_ls = [None]*n_phens
ct_ls = [0]*n_phens
variants_ls = ['']*n_phens
diff_ls = ['']*n_phens
diff_se_ls = ['']*n_phens
diff_pval_ls = ['']*n_phens

for phen_i, phen in enumerate(phens):
    path = wd+f'{phen}.diffgwasloci.tsv.bgz'
    try:
        tb = hl.import_table(path,impute=True)
        tb_ls[phen_i] = tb
        tb1 = tb.filter(tb.diff_pval < threshold).to_pandas()
        ct = len(tb1)
        tb1_ls[phen_i] = tb1
        ct_ls[phen_i] = ct
        variants_ls[phen_i] = ','.join(tb1.variant.values) 
        diff_ls[phen_i] = ','.join(tb1['diff'].values.astype(str)) 
        diff_se_ls[phen_i] = ','.join(tb1['diff_se'].values.astype(str)) 
        diff_pval_ls[phen_i] = ','.join(tb1['diff_pval'].values.astype(str)) 
    except: 
        print(f'\n\nWARNING: Phenotype {phen} ({descs[phen_i]}) is missing results!\n\n')
    print(f'\n\n{phen}, {descs[phen_i]} ({phen_i+1}/{len(phens)}): {ct_ls[phen_i]}\n\n')

#    if phen_i%200==0 or phen_i==len(phens):
df = pd.DataFrame(list(zip(phens,descs,ct_ls,variants_ls,diff_ls,diff_se_ls,diff_pval_ls)),
                  columns=['phenotype','description','sig_loci_ct','sig_loci','diff','diff_se','diff_pval'])

hl.Table.from_pandas(df).export(wd+'ukb31063.diff_gwas_loci.v2.tsv')



# Locally

#tb = hl.import_table('/Users/nbaya/Downloads/ukb31063.diff_gwas_loci.v2.tsv.bgz',force_bgz=True,impute=True)
#
#tb = tb.annotate(sig_loci_ct = hl.int(tb.sig_loci_ct))
#
#tb.export('/Users/nbaya/Downloads/ukb31063.diff_gwas_loci.v2.tsv')


