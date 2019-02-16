# -*- coding: utf-8 -*-
"""
Created on Sun Apr  8 15:21:13 2018

@author: YudongCai

@Email: yudongcai216@gmail.com
"""

import click
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def calTau(exparray: np.array):
    """
    exparray should be log-transformed
    """
    maxval = np.max(exparray)
    n = len(exparray)
    xf = exparray / maxval
    return np.sum(1 - xf) / (n - 1)

def draw(ay, outfile):
    fig, ax = plt.subplots(1, 1, figsize=(8,4))
    bins = ax.hist(ay, bins = 'auto')
    ax.set_xlim(0, 1)
    plt.savefig(outfile)
    plt.close()

    
@click.command()
@click.option('--expfile', help='tab seprated raw fpkm file')
@click.option('--tissuesfile', help='tissues column name use to cal Tau')
@click.option('--indexcol', help='index column name')
@click.option('--outprefix', help='outfile prefix')
def main(expfile, tissuesfile, indexcol, outprefix):
    """
    calculate Tau(tissue specificity).
    input FPKM is raw FPKM
    it will be log-transformed before cal Tau
    out FPKM is log-transformed
    1. Kryuchkova-Mostacci, N. & Robinson-Rechavi, M. A benchmark of gene expression tissue-specificity metrics. Brief. Bioinform. 18, 205–214 (2017).
    """
    tissues = [x.strip() for x in open(tissuesfile)]
    df = pd.read_csv(expfile, sep='\t', usecols=tissues+[indexcol], index_col=indexcol)
    df[df < 1] = 1
    df[tissues] = np.log2(df[tissues].values)
    df['Tau'] = df[tissues].apply(calTau, axis=1)
    df.to_csv(f'{outprefix}.tsv.gz', sep='\t', compression='gzip')
    draw(df['Tau'], f'{outprefix}.tau.pdf')

if __name__ == "__main__":
    main()
