# -*- coding: utf-8 -*-
"""
Created on Mon Jan 15 09:59:34 2018

@author: YudongCai

@Email: yudongcai216@gmail.com
"""

import json
import click
import numpy as np
import pandas as pd
from scipy.stats import norm



def get_cutoff(rvs: list, sperc=0.005):
    """
    输入的rvs是符合正态分布的一组随机变量(random variables)
    sperc, survival percent, 正态分布右侧剩余面积比例, 如果取最大的0.005则设为0.005, 如果取最小的0.005则设为0.995
    """
    mean = np.mean(rvs)
    std = np.std(rvs)
    cutoff = norm.isf(sperc, loc=mean, scale=std)
    return cutoff


def load_data(infile, chr_col, val_col):
    """
    """
    df = pd.read_csv(infile, sep='\t', dtype={chr_col: str, val_col: float})
    return df


@click.command()
@click.option('--infile', help='input tsv data file, first row is header')
@click.option('--chr-col', help='染色体列名', default=False)
@click.option('--sepchrom', is_flag=True, default=False, help='flag, 每条染色体分开算')
@click.option('--sexchrom', help='性染色体名称, defalut is None, 如果使用这个选项, 那么单把性染色体和常染色体分开算, 忽略--sepchrom', multiple=True, default=None)
@click.option('--val-col', help='待计算值列名')
@click.option('--outprefix', help='输出文件前缀')
@click.option('--sperc', help='取top设0.005，取bottom设0.995，默认0.005', default=0.005)
def main(infile, chr_col, sepchrom, sexchrom, val_col, outprefix, sperc):
    """
    apply Inverse survival function to normally distributed variables
    inverse of (1 - cdf)
    """
    print(__doc__)
    print(main.__doc__)
    df = load_data(infile, chr_col, val_col)
    print(f'raw: {df.shape}')
    df.dropna(inplace=True)
    print(f'filter NA, {df.shape}')
    chroms = df[chr_col].unique()
    cutoff = {}
    if sexchrom:
        sexmask = df[chr_col].isin(sexchrom).values
        rvs_auto = df[val_col].values[~sexmask]
        rvs_sex = df[val_col].values[sexmask]
        autocut = get_cutoff(rvs_auto, sperc)
        sexcut = get_cutoff(rvs_sex, sperc)
        for chrom in chroms:
            cutoff[chrom] = autocut if chrom not in sexchrom else sexcut
    elif sepchrom and (not sexchrom):
        for chrom in chroms:
            rvs = df.loc[df[chr_col]==chrom, val_col].values
            cutoff[chrom] = get_cutoff(rvs, sperc)
    elif (not sepchrom) and (not sexchrom):
        cutoff_val = get_cutoff(df[val_col].values, sperc)
        for chrom in chroms:
            cutoff[chrom] = cutoff_val
    else:
        pass
    with open(f'{outprefix}.json', 'w') as f:
        json.dump(cutoff, f, indent=4)
    # 输出显著的记录
    df.loc[df[val_col]>=df[chr_col].map(cutoff), :].to_csv(f'{outprefix}_pass.tsv.gz',
                                                           sep='\t', index=False, compression='gzip')


if __name__ == '__main__':
    main()
