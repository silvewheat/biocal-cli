# -*- coding: utf-8 -*-
"""
Created on Thu Aug  9 21:36:58 2018

@author: YudongCai

@Email: yudongcai216@gmail.com
"""

import click
import numpy as np
import pandas as pd
from scipy.stats import norm



def load_data(infile, val_col):
    """
    """
    df = pd.read_csv(infile, sep='\t', dtype={val_col: float})
    return df


@click.command()
@click.option('--infile', help='input tsv data file, first row is header')
@click.option('--val-col', help='待计算值列名')
@click.option('--outfile', help='输出文件名')
def main(infile, val_col, outfile):
    """
    calculate empirical p value (norm distribution)
    """
    df = load_data(infile, val_col)
    mean = df[val_col].mean()
    std = df[val_col].std()
    df['empirical_Pvalue'] = df[val_col].apply(lambda x: norm.sf(x, loc=mean, scale=std))
    df.to_csv(outfile, sep='\t', index=False)


if __name__ == '__main__':
    main()


