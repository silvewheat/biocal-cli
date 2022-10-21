# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 15:36:21 2021

@author: Yudongcai

@Email: yudong_cai@163.com
"""

import gzip
import typer
import numpy as np
import pandas as pd

def unique(ay):
    return np.unique(ay)[0]


def main(infile: str = typer.Argument(..., help="输入文件名(.anc.tsv.gz)"),
         outfile: str = typer.Argument(..., help='输出文件名(.indseg.tsv.gz)')):
    """合并上一步loter_multiRegions.py的输出结果XX.anc.tsv.gz"""
    cols = pd.read_csv(infile, sep='\t', nrows=0).columns
    hapIDs = cols[3:]
    dtypes = {x: 'category' for x in cols}
    dtypes['pos'] = int # 除了pos是int，别的都是category

    # 先输出个header
    header = '\t'.join(['chrom', 'start', 'end', 'nSNPs', 'sourceIndex', 'hapID']) + '\n'
    with gzip.open(outfile, 'wb') as f:
        f.write(header.encode())

    df = pd.read_csv(infile, sep='\t', dtype=dtypes)

    for hapID in hapIDs:
        cumindex = (df[['chrom', 'region', hapID]] != df[['chrom', 'region', hapID]].shift()).apply(max, axis=1).cumsum()
        fundict = {'chrom': unique, 'pos': [min, max, len], hapID: unique}
        mdf = df[['chrom', 'pos', hapID]].groupby(cumindex).agg(fundict)
        mdf.columns = ['chrom', 'start', 'end', 'nSNPs', 'sourceIndex']
        mdf['hapID'] = hapID
        mdf.to_csv(outfile, sep='\t', index=False, header=None, mode='a', compression='gzip')

if __name__ == '__main__':
    typer.run(main)
    