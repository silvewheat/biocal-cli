import sys
import os
import glob
import numpy as np
import pandas as pd

def unique(ay):
    return np.unique(ay)[0]

def process_df(df):
    cumindex = (df[['hapID', 'sourceIndex']] != df[['hapID', 'sourceIndex']].shift()).apply(max, axis=1).cumsum()
    fundict = {'chrom': unique,
               'pos': [min, max, len],
               'hapID': unique, 
               'group': unique,
               'sourceIndex': unique,
               'nb_bagging': np.mean}
    mdf = df.groupby(cumindex).agg(fundict)
    mdf.columns = ['chrom', 'start', 'end', 'nSNPs', 'hapID', 'group', 'sourceIndex', 'mean_nb_bagging']
    mdf = mdf.loc[mdf['sourceIndex']!=0, :]
    return mdf

refpops = 'bezoar,Capra_falconeri,Capra_caucasica,Capra_nubiana,Capra_ibex,Capra_pyrenaica,Capra_sibirica'
refpops = pd.Series(refpops.split(','))

file = sys.argv[1]
outfile = file[:-7] + '.mdf.tsv.gz'
if not os.path.exists(outfile):
    df = pd.read_csv(file, sep='\t', dtype={'chrom': 'category',
                                            'hapID': 'category',
                                            'group': 'category',
                                            'sourceIndex': 'category'})
    if df.shape[0] > 0:
        df = process_df(df)
        df['source'] = df['sourceIndex'].astype(int).map(refpops)
        df.to_csv(outfile, sep='\t', index=False)