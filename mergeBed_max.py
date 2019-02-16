# -*- coding: utf-8 -*-
"""
Created on Fri Oct  5 20:45:42 2018

@author: YudongCai

@Email: yudongcai216@gmail.com
"""

import os
import click
import numpy as np
import pandas as pd



def load_bed(bedfile):
    """
    四列，没header
    """
    df = pd.read_csv(bedfile, sep='\s+', header=None,
                     names=['chrom', 'start', 'end', 'value'],
                     dtype={'chrom': str,
                           'start': int, 'end': int, 'value': int})
    return df


def agg2df(values):
    last_left = 0
    last_value = 0
    df = []
    for index, value in enumerate(values, 1):
        if value != last_value:
            if last_value > 0:
                df.append([last_left, index-1, last_value])
                last_left = index
                last_value = value
            else:
                last_left = index
                last_value = value
        else:
            last_value = value
    if value == last_value and last_value > 0:
        df.append([last_left, index, last_value])
    df = pd.DataFrame(df, columns=['start', 'end', 'value'])
    return df


def cal_outlier(counts, x=3):
    Q1 = np.percentile(counts, 25)
    Q3 = np.percentile(counts, 75)
    IQR = Q3 - Q1
    upper_boundary = Q3 + IQR * x
    lower_boundary = Q1 - IQR * x
    return upper_boundary, lower_boundary



@click.command()
@click.option('--bedfile', help='四列,没header,第四列是值')
@click.option('--outfile', help='输出文件')
@click.option('--k', help='boxplot排除掉异常高的值, 默认3倍IQR', default=3, type=float)
@click.option('--hardcutoff', help='硬cutoff, 大于这个的去掉, 设了这个会忽略--k', type=float, default=None)
@click.option('--zerobaed', help='flag, 如果是0base,就加上,不加按1base处理', is_flag=True, default=False)
def main(bedfile, outfile, k, hardcutoff, zerobaed):
    try:
        os.remove(outfile)
    except Exception:
        pass
    df = load_bed(bedfile)
    print(f'bedfile loaded {df.shape}')
    # 0base 1base
    if zerobaed:
        df['start'] += 1
        print('0 based')
    else:
        print('1 based')
    # 去异常值
    if not hardcutoff:
        cutoff, _ = cal_outlier(df['value'].values, x=k)
        print(f'length filter cutoff: {cutoff}')
        df = df.loc[(df['end'].values-df['start'].values+1) <= cutoff, :]
        print(f'{df.shape} remained')
    else:
        print(f'cutoff: {hardcutoff}')
        df = df.loc[(df['end'].values-df['start'].values+1) <= hardcutoff, :]
        print(f'{df.shape} remained')
    for chrom in df['chrom'].unique():
        beds = df.loc[df['chrom']==chrom, ['start', 'end', 'value']].values
        length = beds[:, 1].max()
        values = np.zeros(length, dtype=np.int64)
        for start, end, value in beds:
            old = values[start-1: end]
            new = np.full(end-start+1, value)
            values[start-1: end] = np.max([old, new], axis=0)
        odf = agg2df(values)
        odf['chrom'] = chrom
        if zerobaed:
            odf['start'] -= 1
        odf[['chrom', 'start', 'end', 'value']].\
            to_csv(outfile, mode='a', header=False, sep='\t', index=False)


if __name__ == '__main__':
    main()


