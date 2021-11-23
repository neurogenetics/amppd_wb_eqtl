import sys
import pandas as pd
import numpy as np
from numpy import argsort
import dask.dataframe as dd
import pyarrow.parquet as pq
import csv
import time
from gtfparse import read_gtf

import os

import warnings
#warnings.filterwarnings('ignore')
pd.set_option('display.max_columns', None)


cohort = sys.argv[1]
version = sys.argv[2]
tissue = sys.argv[3]
cell = sys.argv[4]

# naming
cohort_version = f'{cohort}.{version}'
cohort_build = f'{cohort}.{version}.{tissue}'

# directories
wrk_dir = f'/labshare/anni/eqtl/tensorqtl_meta/{cohort}'
geno_dir = f'{wrk_dir}/genotypes'
expr_dir = f'{wrk_dir}/expression'
info_dir = f'{wrk_dir}/sample_info'
tensorqtl_dir = f'{wrk_dir}/tensorqtl'
results_dir = f'{wrk_dir}/results'

# input files
bfile_prefix_path = f'{geno_dir}/{cohort_version}.bfile'

def make_esd_files(cohort, version, cohort_build, cell_type):
    print(f'{cell_type}')
    parquet_dir = f'/labshare/anni/eqtl/tensorqtl_meta/{cohort}/tensorqtl'
    cell_files = f'{parquet_dir}/{cohort_build}.{cell_type}.cis_qtl_pairs.chr*.parquet'
    cieqtl_df = dd.read_parquet(cell_files)
    #cieqtl_df['new_gene'] = cieqtl_df['phenotype_id'].str.partition('.')[0]
    probe_ids = list(set(cieqtl_df['phenotype_id']))
    #probe_ids = cieqtl_df['new_gene'].tolist()
    print(f'genes: {len(probe_ids)}')

    #print(cieqtl_df.shape)
    #display(cieqtl_df.head())
    bim_dir = f'/labshare/anni/eqtl/tensorqtl_meta/{cohort}/genotypes'
    genotype_df = dd.read_csv(f'{bim_dir}/{cohort}.{version}.bfile.bim', sep = '\t', header=None)
    genotype_df = genotype_df.rename(columns={0:'chr',1:'variant_id',3:'pos',4:'ref',5:'alt'})

    genotype_df['variant_id'] = genotype_df['variant_id'].str.replace('_b38','')
    genotype_df['variant_id'] = genotype_df['variant_id'].str.replace('_',':')
    #genotype_df.head()

    merge_df = dd.merge(cieqtl_df,genotype_df, on='variant_id')

#make esd files
    print('making esd files...')
    ##subset columns
    esd_df = merge_df[['phenotype_id','chr','variant_id','pos','ref','alt','maf','b_gi','b_gi_se','pval_gi']]
    esd_df = esd_df.drop_duplicates()

    ##convert to numpy array to speed up subsetting
    esd_array = esd_df.compute().to_numpy()
    #new_ids = list(set([item[0] for item in esd_array]))

    print('subsetting genes...')
    for gene in probe_ids:
        ##make subset array per gene
        probe_array = esd_array[np.in1d(esd_array[:,0],gene)]

        ##remove gene id column 
        #probe_array = np.delete(probe_array, 1, 0)
        probe_array = [i[1:] for i in probe_array]
        ##save to textfile
        filename = f'/labshare/anni/eqtl/tensorqtl_meta/{cohort}/osca/{cohort}.{cell_type}.{gene}.esd'
        with open(filename,"w+") as my_csv:
            print('Chr\tSNP\tBp\tA1\tA2\tFreq\tBeta\tse\tp', file=my_csv)
            csvWriter = csv.writer(my_csv,delimiter='\t')
            csvWriter.writerows(probe_array)

    print(f'{cell_type}: done.')

#for cell in cell_types_list:
make_esd_files(cohort, version, cohort, cell)
