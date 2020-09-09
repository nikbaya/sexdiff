#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 14:26:44 2019

Compare results of sex-stratified PRS regression

Does MTAG always outperform

@author: nbaya
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

sexdiff_wd = '/Users/nbaya/Documents/lab/ukbb-sexdiff/'

phen_dict = {
            '50_irnt':['Standing height', 73178],
            '23105_irnt':['Basal metabolic rate', 35705],
            '23106_irnt':['Impedance of the whole body', 73701],
            '2217_irnt':['Age started wearing contact lenses or glasses', 73178],
            '23127_irnt':['Trunk fat percentage', 73178],
            '1100':['Drive faster than motorway speed limit', 73178],
            '1757':['Facial ageing', 35705],
            '6159_3':['Pain type(s) experienced in last month: Neck or shoulder',73178],
            '894':['Duration of moderate activity',73178],
            '1598':['Average weekly spirits intake',35705]
            }

df = pd.read_csv(f'{sexdiff_wd}prs_phen_reg.all_10_phens.n_remove_10000.pruned.tsv',
                 sep='\t')

mtag_def = df[df.gwas_version=='mtag_def'] # adjusted male/female results using MTAG

mtag_rg1 = df[df.gwas_version=='mtag_rg1'] # meta-analyzed male/female results using MTAG

unadjusted = df[df.gwas_version=='unadjusted']

mtag = df[df.gwas_version.str.contains('mtag')]
unadj = df[df.gwas_version=='unadjusted']

merge = mtag.merge(unadj, on=['phen','gwas_sex','sex_tested_on','percentile'],suffixes=['_mtag','_unadjusted'])

merge['multiple_r2_diff'] = merge['multiple_r2_mtag']-merge['multiple_r2_unadjusted']
merge['adjusted_r2_diff'] = merge['adjusted_r2_mtag']-merge['adjusted_r2_unadjusted']

merge['adjusted_r2_diff'].mean()

merge['adjusted_r2_diff'].median()

merge[merge.adjusted_r2_diff<0]

percentile=1
fig, ax = plt.subplots(figsize=(6,4))
for field in ['multiple_r2_diff','adjusted_r2_diff']:
    n, _, _ = ax.hist(merge[merge.percentile==1][field],bins=np.linspace(-0.005,0.05,89),alpha=0.5)
plt.title(f'Distribution of MTAG - unadjusted PRS R^2\n(percentile={percentile})')
plt.legend(['multiple_r2_diff','adjusted_r2_diff'])
ax.plot([0,0],[0, max(n)*1.05],'k--',lw=1)
plt.ylim([0, max(n)*1.05])

for phen, info in phen_dict.items():
    print(f'\n\n... phen: {phen} ({phen_dict[phen][0]}) ...')
    for sex in ['both_sexes','male','female']:
        print(f'\n... gwas sex: {sex} ...')    
        for percentile in [1]: #,0.5,0.1]:
            print(f'\n... percentile: {percentile} ...')    
            print(merge[((merge.phen==phen)&(merge.gwas_sex==sex)&(merge.percentile==percentile))][['phen','sex_tested_on','adjusted_r2_diff']])#'multiple_r2_mtag','multiple_r2_unadjusted']])
