#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 18:01:18 2019

For a given phenotype:
1) Randomly remove 10k males, 10k females from UKB white British with phenotype defined.
2) Run GWAS in males with 10k removed, females with 10k removed for phenotype (using HM3 SNPs).
3) Run GWAS in combined set of males+females (with previously defined 10k males, 10k females removed)
4) Meta-analyze (male, female GWAS?) with MTAG, assuming rg=1
5) Meta-analyze with MTAG without assuming rg=1
6) Create PRS using P+T and each of the 3 GWAS results.
7) Compute difference in PRS-phenotype R2 using PRS from combined GWAS and PRS male GWAS, female GWAS


@author: nbaya
"""

import hail as hl
import numpy as np
from hail.utils.java import Env
import pandas as pd
from scipy import stats
from datetime import datetime as dt
import argparse
from ukbb_pan_ancestry.resources.genotypes import get_ukb_imputed_data # https://github.com/atgu/ukbb_pan_ancestry/tree/master
#import requests
#url = 'https://raw.githubusercontent.com/nikbaya/split/master/gwas.py'
#r = requests.get(url).text
#exec(r)
#gwas=gwas

prs_dir = 'gs://nbaya/sexdiff/prs/'

phen_dict = {
    ## v1 phenotypes
            '50_irnt':['Standing height', 73178],
            '23105_irnt':['Basal metabolic rate', 35705],
            '23106_irnt':['Impedance of the whole body', 73701],
            '2217_irnt':['Age started wearing contact lenses or glasses', 73178],
            '23127_irnt':['Trunk fat percentage', 73178],
            '1100':['Drive faster than motorway speed limit', 73178],
            '1757':['Facial ageing', 35705],
            '6159_3':['Pain type(s) experienced in last month: Neck or shoulder',73178],
            '894':['Duration of moderate activity',73178],
            '1598':['Average weekly spirits intake',35705],
    ## raw v1 phenotypes
            '50_raw':['Standing height (raw)', 26308],
            '23105_raw':['Basal metabolic rate (raw)', 73178],
            '23106_raw':['Impedance of the whole body (raw)', 73178],
            '2217_raw':['Age started wearing contact lenses or glasses (raw)', 73178],
            '23127_raw':['Trunk fat percentage (raw)', 73178],
    ## v2 phenotypes
              'K11_HERNIA': ['Hernia', 73178],
              '4079_irnt': ['Diastolic blood pressure, automated reading', 73178],
              '4080_irnt': ['Systolic blood pressure, automated reading', 73178],
        }

sexes = ['female','male','both_sexes']

n_remove_per_sex = 10e3 # number of samples to remove for each sex
prune=True # whether to LD prune variants when calculating PRS, default: True
percentiles = [1, 0.5, 0.1] # p-value percentiles, default: [1, 0.5, 0.1]
use_sex_spec_irnt = True # whether to use sex-specific phenotype information when comparing against PRS


def get_mt(overwrite=False):
    mt_path = 'gs://ukbb-temp-30day/nbaya/ukb31063.hm3_variants.gwas_samples.mt'
    if hl.hadoop_is_file(f'{mt_path}/_SUCCESS'):
        mt = hl.read_matrix_table(mt_path)
    else:
        tmp_mt_path = mt_path.replace('nbaya/','nbaya/tmp-')
        def _write_tmp_mt():
            ## part 1: about an 1 hr with 30 workers (possibly starting with 10 then increasing to 30 if progress stalls) 
            variants = hl.import_table('gs://nbaya/hapmap3_variants.tsv.gz', force=True) # download here: https://github.com/nikbaya/split/blob/master/hapmap3_variants.tsv.gz
            variants = variants.key_by(**hl.parse_variant(variants.v))
            mt = get_ukb_imputed_data('all', variant_list = variants, entry_fields = ('dosage', )) # 'all' = autosomes only
            # print(mt.count()) # (1089172, 487409)        
            # mt = mt.checkpoint(mt_path.replace('nbaya/','nbaya/tmp-'), overwrite=overwrite)
            mt.write(tmp_mt_path, _read_if_exists=True)
        
        def _repartition():
            ## part 2: 5 min with 100 preemptibles
            mt = hl.read_matrix_table(tmp_mt_path)
            mt = mt.repartition(1000)
            withdrawn = hl.read_table('gs://ukb31063/ukb31063.withdrawn_samples.ht')
            mt = mt.anti_join_cols(withdrawn)
            # print(mt.count()) # (1089172, 487409)
            covs = hl.read_table('gs://ukb31063/ukb31063.neale_gwas_covariates.both_sexes.ht')
            mt = mt.annotate_cols(**covs[mt.s])
            mt = mt.filter_cols(hl.is_defined(mt.PC1))
            # print(mt.count())  # (1089172, 361144)
            return mt.checkpoint(mt_path, overwrite=overwrite)
        _write_tmp_mt()
        mt = _repartition()
        
    return mt

def remove_n_individuals(mt, n_remove_per_sex, phen, phen_tb_dict, sexes = ['females','males'], seed=None):
    r'''
    Removes n_remove_per_sex individuals from each specified sex (default is to remove from both
    females and males, i.e. sexes='fm').
    Saves tables with phenotype data for each sex and the sexes combined.
    '''
    assert 'female' in sexes or 'male' in sexes, "sexes must have 'female' or 'male'"
    n_remove_per_sex = int(n_remove_per_sex)
    print(f'\n... Removing {n_remove_per_sex} per sex from UKB phenotype ({phen}) ...\n')
    seed = seed if seed is not None else int(str(Env.next_seed())[:5])
    print(f'\n... seed: {seed} ...\n')
    mt_cols = mt.cols()
#    initial_ct = mt_cols.count()
#    print(f'\n\n... Initial combined mt count ({phen}): {initial_ct} ...\n')
    tb_sex_ls = [None, None]
    for idx, sex in enumerate(sexes):
        tb_sex_path = prs_dir+f'{phen}.{sex}.n_remove_{n_remove_per_sex}.seed_{seed}.ht'
        if hl.hadoop_is_file(f'{tb_sex_path}/_SUCCESS'):
            print(f'\n... {phen} {sex} table already written! ...\n')
        else:
            print(f'\n... Starting to write {phen} table for {sex} ...\n')
            mt_cols_sex = annotate_phen(tb=mt_cols, phen=phen, sex=sex, phen_tb_dict=phen_tb_dict)
            n_sex = mt_cols_sex.count()
            print(f'\n\n... Initial {sex} mt count ({phen}): {n_sex} ...\n')
            tb1 = mt_cols_sex.add_index('idx_tmp')
            tb2 = tb1.key_by('idx_tmp')
            remove = [1]*(n_remove_per_sex)+[0]*(n_sex-n_remove_per_sex)
            randstate = np.random.RandomState(int(seed)) 
            randstate.shuffle(remove)
            tb3 = tb2.annotate(remove = hl.literal(remove)[hl.int32(tb2.idx_tmp)])
            tb4 = tb3.filter(tb3.remove == 1, keep=True) #keep samples we wish to discard from original mt
            tb5 = tb4.key_by('s')
            tb5.select('phen').write(tb_sex_path,overwrite=True) # write out table with samples of single sex we wish to discard from original mt
        tb_sex_ls[idx] = hl.read_table(tb_sex_path)
    if len(set(sexes).union(['female','male'])) == 2: #if sexes has both sexes
        tb_both_path = prs_dir+f'{phen}.both_sexes.n_remove_{n_remove_per_sex}.seed_{seed}.ht'
        if hl.hadoop_is_file(f'{tb_both_path}/_SUCCESS'):
            print(f'\n... {phen} both_sexes table already written! ...\n')
        else:
            print(f'\n... Starting to write {phen} table for both_sexes ...\n')
            mt_cols_both = annotate_phen(tb=mt_cols, phen=phen, sex='both_sexes',
                                         phen_tb_dict=phen_tb_dict)
            n_both = mt_cols_both.count()
            print(f'\n\n... Initial both_sexes mt count ({phen}): {n_both} ...\n')
            tb_both1 = mt_cols_both.anti_join(tb_sex_ls[0]) # remove individuals in sex1 holdout set
            tb_both2 = tb_both1.anti_join(tb_sex_ls[1]) # remove individuals in sex2 holdout set
            tb_both = mt_cols_both.anti_join(tb_both2) # remove individuals we wish to KEEP from mt_cols_both, leaving only those we are discarding
            tb_both.select('phen').write(tb_both_path,overwrite=True) # write out table containing all m/f individuals we wish to discard from original mt
        tb_both = hl.read_table(tb_both_path)
        ht_both = annotate_phen(tb=mt_cols, phen=phen, sex='both_sexes', phen_tb_dict=phen_tb_dict)
        mt_both = mt.annotate_cols(phen = ht_both[mt.s].phen)
        mt_both = mt_both.filter_cols(hl.is_defined(mt_both.phen))
        mt_both = mt_both.anti_join_cols(tb_both)
        mt_both_ct = mt_both.count_cols()
        print(f'\n\n... Final both_sexes mt count ({phen}): {mt_both_ct} ...\n')
    else: 
        mt_both=None
    ht_f = annotate_phen(tb=mt_cols, phen=phen, sex='female', phen_tb_dict=phen_tb_dict)
    mt_f = mt.annotate_cols(phen = ht_f[mt.s].phen)
    mt_f = mt_f.filter_cols(hl.is_defined(mt_f.phen))
    mt_f = mt_f.anti_join_cols(tb_sex_ls[sexes.index('female')]) if mt_f is not None else None
    ht_m= annotate_phen(tb=mt_cols, phen=phen, sex='male', phen_tb_dict=phen_tb_dict)
    mt_m = mt.annotate_cols(phen = ht_m[mt.s].phen)
    mt_m = mt_m.filter_cols(hl.is_defined(mt_m.phen))
    mt_m = mt_m.anti_join_cols(tb_sex_ls[sexes.index('male')]) if mt_m is not None else None
    mt_sex_ls = [mt_f, mt_m]
    for idx, sex in enumerate(sexes):
        mt_sex_ct = mt_sex_ls[idx].count_cols()
        print(f'\n\n... Final {sex} mt count ({phen}): {mt_sex_ct} ...\n')

    return mt_both, mt_f, mt_m, seed
    
def annotate_phen(tb, phen, sex, phen_tb_dict, filter_to_phen=True):
    r'''
    Annotates `tb` with phenotype `phen` and filters to individuals with 
    phenotype defined. Uses sex-specific IRNT phenotypes.
    sex options: female, male, both_sexes
    '''
    print(f'\n... Reading UKB phenotype "{phen_dict[phen][0]}" for {sex} (code: {phen}) ...')
        
    phen_tb0 = phen_tb_dict[sex]
    phen_tb = phen_tb0.select(phen).rename({phen:'phen'})

    if type(tb)==hl.table.Table:
        annotate_fn = hl.Table.annotate
        filter_fn = hl.Table.filter
    elif type(tb)==hl.matrixtable.MatrixTable:
        annotate_fn = hl.MatrixTable.annotate_cols
        filter_fn = hl.MatrixTable.filter_cols

    tb0 = annotate_fn(self=tb, phen_str = hl.str(phen_tb[tb.s]['phen']).replace('\"',''))
    
    if filter_to_phen: # filter to individuals with phenotype data defined
        tb0 = filter_fn(self=tb0, expr=tb0.phen_str == '', keep=False)
    
    if phen_tb.phen.dtype == hl.dtype('bool'):
        tb1 = annotate_fn(self=tb0, phen = hl.bool(tb0.phen_str)).drop('phen_str')
    else:
        tb1 = annotate_fn(self=tb0, phen = hl.float64(tb0.phen_str)).drop('phen_str')
    
    return tb1

def gwas(mt, x, y, cov_list=[], with_intercept=True, pass_through=[], path_to_save=None, 
         normalize_x=False, is_std_cov_list=True):
    r'''Runs GWAS'''
    
    mt = mt._annotate_all(col_exprs={'__y':y},
                           entry_exprs={'__x':x})
    
    print('\n... Calculating allele frequency ...')
    mt_freq_rows = mt.annotate_rows(freq = hl.agg.mean(mt.dosage)/2).rows() #frequency of alternate allele
    mt_freq_rows = mt_freq_rows.key_by('rsid')
    
    if normalize_x:
        mt = mt.annotate_rows(__gt_stats = hl.agg.stats(mt.__x))
        mt = mt.annotate_entries(__x= (mt.__x-mt.__gt_stats.mean)/mt.__gt_stats.stdev) 
        mt = mt.drop('__gt_stats')
    
    if is_std_cov_list:
        cov_list = ['isFemale','age','age_squared','age_isFemale',
                    'age_squared_isFemale']+['PC{:}'.format(i) for i in range(1, 21)]
        
    if str in list(map(lambda x: type(x),cov_list)):
        cov_list = list(map(lambda x: mt[x] if type(x) is str else x,cov_list))
        
    cov_list = ([1] if with_intercept else [])+cov_list
    
    print(f'pass through: {pass_through}')

    gwas_ht = hl.linear_regression_rows(y=mt.__y,
                                        x=mt.__x,
                                        covariates=cov_list,
                                        pass_through = ['rsid']+pass_through)
    
    gwas_ht = gwas_ht.annotate_globals(with_intercept = with_intercept)
        
    gwas_ht = gwas_ht.key_by('rsid')
    
    ss_template = hl.read_table('gs://nbaya/rg_sex/hm3.sumstats_template.ht') # sumstats template as a hail table
    ss_template = ss_template.key_by('SNP')
        
    ss = ss_template.annotate(chr = gwas_ht[ss_template.SNP].locus.contig,
                              bpos = gwas_ht[ss_template.SNP].locus.position,
                              freq = mt_freq_rows[ss_template.SNP].freq,
                              beta = gwas_ht[ss_template.SNP].beta,
                              z = gwas_ht[ss_template.SNP].t_stat,
                              pval = gwas_ht[ss_template.SNP].p_value,
                              n = gwas_ht[ss_template.SNP].n)
    ss = ss.drop('N')
    ss = ss.rename({'SNP':'snpid',
                    'A1':'a1',
                    'A2':'a2'})
    
    print(ss.describe())
    
    
    if path_to_save is not None:
        ss.export(path_to_save)
        
    return ss

def get_freq(mt, phen, sex, n_remove, seed):
    r'''
    Get allele frequencies and other SNP information (needed to fix previously 
    created sumstats files)
    '''
    
    print('... Calculating allele frequency ...')
    mt = mt.annotate_rows(freq = hl.agg.mean(mt.dosage)/2) #frequency of alternate allele
    mt_rows = mt.rows()
    mt_rows = mt_rows.key_by('rsid')
    mt_rows = mt_rows.annotate(chr = mt_rows.locus.contig,
                               bpos = mt_rows.locus.position)

    
    ss = hl.import_table(prs_dir+f'{phen}.gwas.{sex}.n_remove_{n_remove_per_sex}.seed_{seed}.old.tsv.bgz',
                         impute=True,
                         key='SNP')
    
    ss = ss.annotate(chr = mt_rows[ss.SNP].chr,
                     bpos = mt_rows[ss.SNP].bpos,
                     freq = mt_rows[ss.SNP].freq,
                     z = ((-1)*(ss.beta<0)*hl.abs(hl.qnorm(ss.p_value/2))+
                          (ss.beta>0)*hl.abs(hl.qnorm(ss.p_value/2)))
                     )
                     

    if 'N' in ss.row:
        if 'n' not in ss.row:
            ss = ss.annotate(n = ss.N)
        ss = ss.drop('N')
        
    ss = ss.rename({'SNP':'snpid',
                    'A1':'a1',
                    'A2':'a2',
                    'p_value':'pval'})
    
    ss = ss.key_by()
    ss = ss.select('snpid','chr','bpos','a1','a2','freq','beta','z','pval','n')
    ss = ss.key_by('snpid')
    
    ss.export(prs_dir+f'{phen}.gwas.{sex}.n_remove_{n_remove_per_sex}.seed_{seed}.tsv.bgz')
    
def get_freq_alt(mt, phen, sex, n_remove, seed):
    
    ss = hl.import_table(prs_dir+f'{phen}.gwas.{sex}.n_remove_{n_remove_per_sex}.seed_{seed}.tsv.bgz',
                         impute=True,
                         key='snpid')
    ss = ss.annotate(z = ((-1)*(ss.beta<0)*hl.abs(hl.qnorm(ss.pval/2))+
                          (ss.beta>0)*hl.abs(hl.qnorm(ss.pval/2)))
                     )
    
    if 'N' in ss.row:
        if 'n' not in ss.row:
            ss = ss.annotate(n = ss.N)
        ss = ss.drop('N')
    
    ss = ss.key_by()
    ss = ss.select('snpid','chr','bpos','a1','a2','freq','beta','z','pval','n')
    ss = ss.key_by('snpid')
    
    ss.export(prs_dir+f'{phen}.gwas.{sex}.n_remove_{n_remove_per_sex}.seed_{seed}.tsv.bgz')
    
def prs(mt, phen, sex, n_remove, prune, percentiles, seed, count=True):
    r'''
    Calculate PRS using betas from both sexes and sex-stratified GWAS, as well
    as MTAG meta-analyzed betas. PRS are always calculated on both sexes, 
    regardless of the sex the GWAS was run on.
    P-value thresholds are determined by percentile.
    Set `count`=True if running this for the first time, to be sure that 
    numbers make sense. To speed up, set count_rows=False
    '''
    assert sex in ['both_sexes','female','male'], f'WARNING: sex={sex} not allowed. sex must be one of the following: both_sexes, female, male'
    
    # "def" uses the MTAG results created by using the default settings
    # "rg1" uses the MTAG results created by using the --perfect-gencov flag
    gwas_versions = ['unadjusted',f'mtag_{"rg1" if sex=="both_sexes" else "def"}'] 
    
    mt = mt
        
    for gwas_version in gwas_versions:
        print(f'\n... Calculating PRS for "{phen_dict[phen][0]}" {sex} {gwas_version} ...\n')
        gwas_version_suffix = "" if gwas_version=='unadjusted' else '.'+gwas_version
        gwas_path = (prs_dir+f'{phen}.gwas.{sex}.n_remove_{n_remove_per_sex}.seed_{seed}{gwas_version_suffix}.tsv.{"b" if gwas_version=="unadjusted" else ""}gz')

        ss=hl.import_table(gwas_path,
                           impute=True,
                           key='snpid' if gwas_version=='unadjusted' else 'SNP',
                           force=True)
        
        if prune:
            print('\n... Pruning SNPs ...\n')
            # define the set of SNPs
            pruned_snps_path = 'gs://nbaya/risk_gradients/ukb_imp_v3_pruned.bim' #from Robert Maier (pruning threshold r2=0.2, random 10k UKB sample), download here: https://github.com/nikbaya/split/blob/master/ukb_imp_v3_pruned.bim.gz
            variants = hl.import_table(pruned_snps_path, delimiter='\t', no_header=True, impute=True)
            print(f'\n... Pruning to variants in {pruned_snps_path} ...\n')
            variants = variants.rename(
                {'f0': 'chr', 'f1': 'rsid', 'f3': 'pos'}).key_by('rsid')
#            mt = mt.key_rows_by('rsid')
            # filter to variants defined in variants table
            ss = ss.filter(hl.is_defined(variants[ss['snpid' if gwas_version=='unadjusted' else 'SNP']]))
            if count:
                ct_rows = ss.count()
                print(f'\n\n... SNP count after pruning filter: {ct_rows} ...\n')
        else:
            print(f'\n... Not pruning because prune={prune} ...\n')
            
        for percentile in percentiles:

            ## use path without threshold to check if PRS was already run because it doesn't require calculating the pval threshold
            prs_path_without_threshold = (prs_dir+f'prs.{phen}.{sex}.n_remove_{int(n_remove_per_sex)}.seed_{seed}.{gwas_version}.{"" if prune else "not_"}pruned.pval_thresh_*.perc_{percentile}.tsv')
            # prs_path_without_threshold = (prs_dir+f'prs.{phen}.{sex}.n_remove_{int(n_remove_per_sex)}.seed_{seed}.{gwas_version}.{"" if prune else "not_"}pruned.pval_thresh_*.perc_{percentile}.opposite_sex.tsv')

            if len(hl.hadoop_ls(prs_path_without_threshold))>0:
                print(f'\n\n... Calculation of PRS for "{phen_dict[phen][0]}" {sex} {gwas_version} for percentile {percentile} already completed! ...\n')
            else:
                start = dt.now()
                
                if percentile != 1:
                    threshold = ss.aggregate(hl.agg.approx_quantiles(ss[('' if gwas_version=='unadjusted' else 'mtag_')+'pval'],percentile))
                    ss = ss.filter(ss[('' if gwas_version=='unadjusted' else 'mtag_')+'pval']<=threshold)
                else:    
                    threshold=1
                threshold_str = '{:.4e}'.format(threshold)
                prs_path = prs_dir+f'prs.{phen}.{sex}.n_remove_{int(n_remove_per_sex)}.seed_{seed}.{gwas_version}.{"" if prune else "not_"}pruned.pval_thresh_{threshold_str}.perc_{percentile}.tsv'
                    
                print(f'\n\n... Using p-value threshold of {threshold} for percentile {percentile} ...\n')                    
                mt = mt.annotate_rows(beta = ss[mt.rsid]['beta' if gwas_version=='unadjusted' else 'mtag_beta'])
    
                if count:
                    if percentile != 1:
                        threshold_ct = mt.filter_rows(hl.is_defined(mt.beta)).count_rows()
                    else:
                        threshold_ct = ct_rows
                        
                    print(f'\n\n... Variants remaining after thresholding filter: {threshold_ct} ...\n')
                
                mt = mt.annotate_cols(prs = hl.agg.sum(mt.dosage*mt.beta))
                
                if count:
                    mt_cols_ct = mt.filter_cols(hl.is_defined(mt.prs)).count_cols()
                    
                    print(f'\n\n... Samples with PRS: {mt_cols_ct} ...\n')

                mt.cols().describe()
                
                mt.cols().select('prs').export(prs_path)
                
                elapsed = dt.now()-start
                
                print(f'\n\n... Completed calculation of PRS for "{phen_dict[phen][0]}" {sex} {gwas_version} ...')
                print(f'\n... Elapsed time: {round(elapsed.seconds/60, 2)} min ...\n')


def get_test_mt(mt_all, phen, sex, seed):
    assert sex in ['both_sexes','female','male'], f'WARNING: sex={sex} not allowed. sex must be one of the following: both_sexes, female, male'
    path = prs_dir+f'{phen}.{sex}.n_remove_{int(n_remove_per_sex)}.seed_{seed}.ht'
    ht_sex = hl.read_table(path)
    return mt_all.semi_join_cols(ht_sex)


def prs_phen_reg(test_mt, phen, sex, n_remove, prune, percentiles, seed, 
                 use_sex_spec_irnt=False, overwrite=False):
    
    if use_sex_spec_irnt and 'irnt' not in phen:
        print(f'NOTE: Setting use_sex_spec_irnt=False because phen {phen} is not IRNT')

    test_ht = test_mt.cols()
    
    reg_path = prs_dir+f'prs_phen_reg.{phen}.{sex}.n_remove_{int(n_remove_per_sex)}.seed_{seed}.{"" if prune else "not_"}pruned{".sexspecirnt" if use_sex_spec_irnt else ""}.tsv'

    if hl.hadoop_is_file(reg_path) and not overwrite:
        print(f'... Phen ~ PRS + covariates regression already complete for all gwas versions & percentiles of {phen} {sex} {"sex_spec_irnt" if use_sex_spec_irnt else ""}! ...')
    else:
        
        row_struct_ls = []
        
        gwas_versions = ['unadjusted', f'mtag_{"def" if sex!="both_sexes" else "rg1"}']
        
        for gwas_version in gwas_versions:
            for percentile in percentiles:
                prs_path_without_threshold = prs_dir+f'prs.{phen}.{sex}.n_remove_{int(n_remove_per_sex)}.seed_{seed}.{gwas_version}.{"" if prune else "not_"}pruned*.perc_{percentile}.tsv'
                print(prs_path_without_threshold)
                ls = hl.hadoop_ls(prs_path_without_threshold)
                print('WARNING: More than one file matches {prs_path_without_threshold}' if len(ls)>1 else '')
                print('\n'.join([x['path'] for x in ls]))
                prs_path = ls[0]['path'] # default to using the first path if more than one path exists for a given p-value percentile
                pval_thresh = prs_path.split('pval_thresh_')[1].split('.perc')[0] #previously used for the both_sexes prs
                print(f'... {phen} {sex} {gwas_version} percentile={percentile} ...')
                print(f'... using {prs_path} ...')
                print(f'... pval threshold: {pval_thresh} ...')
                prs_ht = hl.import_table(prs_path,impute=True, key='s',types={'s': hl.tstr})
                test_ht = test_ht.annotate(prs = prs_ht[test_ht.s].prs)
                
                cov_list = ['prs','age','age_squared']+['PC{:}'.format(i) for i in range(1, 21)]
                for isFemale in [0,1]: # test in males, then females
                    test_ht_sex = test_ht.filter(test_ht.isFemale==isFemale)
                    reg = test_ht_sex.aggregate(hl.agg.linreg(y=test_ht_sex.phen, 
                                                              x=[1]+list(map(lambda x: test_ht_sex[x] if type(x) is str else x,cov_list))))
                    print(f'\n\n... {phen} {sex} {gwas_version} percentile={percentile} '+
                          f'applied to {"fe" if isFemale else ""}males {"using sex-spec irnt" if use_sex_spec_irnt else ""} ...\n'+
                          f'\n... multiple R^2: {reg.multiple_r_squared} ...'+
                          f'\n... pval for multiple R^2: {reg.multiple_p_value} ...'+
                          f'\n... adjusted R^2: {reg.adjusted_r_squared} ...')
                    row_struct_ls.append({'phen':phen, 
                                          'gwas_sex':sex, 
                                          'gwas_version':gwas_version,
                                          'sex_spec_irnt':str(use_sex_spec_irnt),
                                          'percentile':str(percentile),
                                          'pval_threshold':pval_thresh,
                                          'sex_tested_on': f'{"fe" if isFemale else ""}males',
                                          'multiple_r2':str(reg.multiple_r_squared),
                                          'multiple_r2_pval':str(reg.multiple_p_value),
                                          'adjusted_r2':str(reg.adjusted_r_squared)})

        ht = hl.Table.parallelize(hl.literal(row_struct_ls,'array<struct{phen: str, gwas_sex: str, gwas_version: str, sex_spec_irnt: str, percentile: str, pval_threshold: str, sex_tested_on: str, multiple_r2: str, multiple_r2_pval: str, adjusted_r2: str}>'))
        ht = ht.annotate(percentile = hl.float(ht.percentile),
                         pval_threshold = hl.float(ht.pval_threshold),
                         multiple_r2 = hl.float(ht.multiple_r2),
                         multiple_r2_pval = hl.float(ht.multiple_r2_pval),
                         adjusted_r2 = hl.float(ht.adjusted_r2))
        ht.show(12)    

        print(f'\n\n... Writing regression results to {reg_path} (overwrite={overwrite})...')
        ht.export(reg_path)

def create_full_results_file(prune, overwrite=False):
    r'''
    Concatenates PRS-phentype regression results into a single table.
    '''
    reg_path_regex = prs_dir+f'prs_phen_reg.*.*.n_remove_{int(n_remove_per_sex)}.seed_*.{"" if prune else "not_"}pruned*.tsv'
    ls = hl.hadoop_ls(reg_path_regex)
    reg_paths = sorted([f['path'] for f in ls])
    df_list = []
    for reg_path in reg_paths:
        with hl.hadoop_open(reg_path) as f:
            df_list.append(pd.read_csv(f, sep='\t'))
    df = pd.concat(df_list, sort=False)
    df.insert(1, 'phen_desc', df.phen.astype(str).apply(lambda x: phen_dict[x][0])) # add phenotype description to dataframe

    all_reg_results_path = prs_dir+f'prs_phen_reg.all_phens.n_remove_{int(n_remove_per_sex)}.{"" if prune else "not_"}pruned.tsv'
    if hl.hadoop_is_file(all_reg_results_path) and not overwrite:
        print('\n... Full PRS-phen regression results already written! ...')
        print(all_reg_results_path)
    else:
        print('\n... Writing PRS-phen regression results ...')
        print(all_reg_results_path)
        with hl.hadoop_open(all_reg_results_path, 'w') as f:
            df.to_csv(f, sep='\t', index=False)
            
def main(args):
    hl.init(log='/tmp/hail.log')

    run_gwas = args.run_gwas
    calc_prs = args.calc_prs
    phen_prs_reg = args.phen_prs_reg
    
    mt0 = get_mt()
    
    if run_gwas or phen_prs_reg:
        phen_tb_dict = {'female':None,'male':None,'both_sexes':None}
        for sex in sexes:
            phen_tb = hl.import_table(f'gs://ukb31063/ukb31063.PHESANT_January_2019.{sex}.tsv.bgz', 
                                      missing='',impute=True,types={'s': hl.tstr}, key='s')
            if 'K11_HERNIA' in phen_dict:
                finngen_tb = hl.import_table(f'gs://ukb31063/ukb31063.FINNGEN_phenotypes.{sex}.tsv.bgz', 
                                             missing='',impute=True,types={'s': hl.tstr}, key='s')
                phen_tb = phen_tb.join(finngen_tb, how='outer')
            phen_tb_dict[sex] =  phen_tb

    for phen, phen_desc in phen_dict.items():

        if run_gwas:
            mt_both, mt_f, mt_m, seed = remove_n_individuals(mt=mt0, n_remove_per_sex=n_remove_per_sex, 
                                                              phen=phen, phen_tb_dict=phen_tb_dict, 
                                                              sexes = ['female','male'], seed=phen_dict[phen][1])
            for mt_tmp, sex in [(mt_f,'female'), (mt_m,'male'), (mt_both,'both_sexes')]:            
                if sex not in sexes:
                    continue
                gwas_path = prs_dir+f'{phen}.gwas.{sex}.n_remove_{n_remove_per_sex}.seed_{seed}.tsv.bgz'
                if hl.hadoop_is_file(gwas_path):
                    print(f'\n... {sex} GWAS on "{phen_dict[phen][0]}" (code: {phen}) already complete! ...\n')
                else:
                    print(f'\n... Running {sex} GWAS on "{phen_dict[phen][0]}" (code: {phen}) ...\n')
                    gwas(mt=mt_tmp, 
                          x=mt_tmp.dosage, 
                          y=mt_tmp['phen'], 
                          path_to_save=gwas_path,
                          is_std_cov_list=True)       
                    
        ## Run MTAG using a cloud VM
        # =============================================================================
        # conda create -n mtag -y -q  python=2.7 scipy numpy pandas bitarray joblib
        # source activate mtag
        # 
        # git clone https://github.com/omeed-maghzian/mtag.git
        # cd mtag
        # gsutil cp gs://nbaya/sexdiff/prs/run_mtag.sh ~/ # download from https://github.com/nikbaya/sexdiff/blob/master/run_mtag.sh if necessary
        # chmod +x ~/run_mtag.sh
        # mkdir ~/data/
        # =============================================================================

        ## Then run run_mtag.sh, passing the phenotype and seed as appropriate

        if calc_prs:
            mt_sex0 = get_test_mt(mt_all=mt0, phen=phen, sex='both_sexes', seed=phen_dict[phen][1])
            for sex in sexes:   
                prs(mt=mt_sex0, phen=phen, sex=sex, n_remove=n_remove_per_sex,
                    prune=prune, percentiles=percentiles, seed=phen_dict[phen][1], count=True)
            
        if phen_prs_reg:
            mt_sex0 = get_test_mt(mt_all=mt0, phen=phen, sex='both_sexes', seed=phen_dict[phen][1])
            for sex in sexes:
                
                if use_sex_spec_irnt and sex!='both_sexes' and 'irnt' in phen:
                    mt_sex1 = annotate_phen(tb=mt_sex0, phen=phen, sex='female', 
                                            phen_tb_dict=phen_tb_dict, filter_to_phen=False)
                    mt_sex2 = mt_sex1.rename({'phen':'phen_female'})
                    mt_sex3 = annotate_phen(tb=mt_sex2, phen=phen, sex='male', 
                                            phen_tb_dict=phen_tb_dict, filter_to_phen=False)
                    mt_sex4 = mt_sex3.rename({'phen':'phen_male'})
                    mt_sex = mt_sex4.annotate_cols(phen = hl.or_else(mt_sex4.phen_female, mt_sex4.phen_male))
                    
                else:
                    mt_sex = annotate_phen(tb=mt_sex0, phen=phen, sex='both_sexes', phen_tb_dict=phen_tb_dict) #use the both sexes irnt
                    
                prs_phen_reg(test_mt=mt_sex, phen=phen, sex=sex, n_remove=n_remove_per_sex, 
                              prune=prune, percentiles=percentiles, seed=phen_dict[phen][1],
                              use_sex_spec_irnt=use_sex_spec_irnt)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--run_gwas',action='store_true',default=False,help="whether to run GWAS")
    parser.add_argument('--calc_prs',action='store_true',default=False,help="whether to calculate PRS")
    parser.add_argument('--phen_prs_reg',action='store_true',default=False,help="whether to run phenotype ~ PRS + covs regression")
    parser.add_argument('--create_full_results_file',action='store_true',default=False,help="whether to create a full results file")
    parser.add_argument('--overwrite',action='store_true',default=False,help="whether to overwrite files")
    args = parser.parse_args()
    
    if args.create_full_results_file:
        create_full_results_file(prune=prune, 
                                 overwrite=args.overwrite)
    else:
        main(args)
    


