# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 11:49:46 2020

@author: YudongCai

@Email: yudongcai216@gmail.com
"""

import click
import numpy as np
import pandas as pd


@click.command()
@click.option('--infile', help='rfmix输出的.msp.tsv文件,分染色体的')
@click.option('--outprefix', help='输出文件前缀')
def main(infile, outprefix):
    """
    输入文件得是分染色体的
    """
    with open(infile) as f:
        line1 = f.readline()
        line2 = f.readline()[1:] # 去掉开头井号
    groupidex2name = {}
    for item in line1.split(':')[1].strip().split('\t'):
        group_name, group_index = item.split('=')
        group_index = int(group_index)
        groupidex2name[group_index] = group_name
    header = line2.strip().split('\t')
    header[5] = 'n_snps' # 强迫症不喜欢名字中间夹空格的
    hapIDs = header[6:]
    df = pd.read_csv(infile, sep='\t', skiprows=2, header=None, names=header)
    df = df.melt(id_vars=header[:6], value_vars=hapIDs, var_name='hapID', value_name='sourceID')
    # msp文件的坐标是连续无gap的，所以：
    cumindex = (df[['hapID', 'sourceID']] != df[['hapID', 'sourceID']].shift()).apply(max, axis=1).cumsum()
    fundict = {'chm': max, 'spos': min, 'epos': max,
               'n_snps': sum, 'hapID': max, 'sourceID': max}
    mdf = df.groupby(cumindex).agg(fundict)
    mdf['length'] = mdf['epos'] - mdf['spos'] + 1
    mdf['sourceName'] = mdf['sourceID'].map(groupidex2name)
    mdf[['chm', 'spos', 'epos', 'length', 'n_snps', 'hapID', 'sourceName']]\
        .to_csv(f'{outprefix}.indSeg.tsv.gz', compression='gzip', index=False, sep='\t')
    ddf = mdf[['hapID', 'sourceName', 'length']].groupby(['hapID', 'sourceName']).describe().reset_index()
    ddf.columns = ['hapID', 'sourceName', 'count', 'mean', 'std', 'min', '25%', '50%', '75%', 'max']
    ddf.to_csv(f'{outprefix}.indstat.tsv.gz', compression='gzip', index=False, sep='\t')


if __name__ == '__main__':
    main()

