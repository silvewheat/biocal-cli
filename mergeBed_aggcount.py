# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 16:46:09 2019

@author: YudongCai

@Email: yudongcai216@gmail.com
"""

import click
import numpy as np
import pandas as pd



@click.command()
@click.option('--infile', help='输入四列bed文件,第四列是计数, 无header')
@click.option('--outfile', help='输出文件')
def main(infile, outfile):
    """
    \b
    合并 chrom start end count 类似文件
    保证chrom和count一致情况下合并, 区间重合的累加count
    坐标系统假设为1-based
    如：
    infile:
    1 101 500 1
    1 201 300 2
    1 401 600 1
    outfile:
    1 101 200 1
    1 201 300 3
    1 301 400 1
    1 401 500 2
    1 501 600 1
__________________
    | |  2
______________
    |||  3
    |||  3
______________
    ||   2
______________
    |    1
__________________
    """
    with open(outfile, 'w') as f:
        df = pd.read_csv(infile, sep='\s+', header=None, names=['chrom', 'start', 'end', 'count'])
        chroms = df['chrom'].unique()
        for chrom in chroms:
            print(chrom)
            tdf = df.loc[df['chrom']==chrom, :]
            length = tdf['end'].max()
            pseudo_chrom = np.zeros(length)
            for chrom, start, end, count in tdf.values:
                pseudo_chrom[start-1: end] += count
            last_count = pseudo_chrom[0]
            if last_count != 0:
                region_start = 0
            for index, count in enumerate(pseudo_chrom, 0):
                if count == last_count:
                    pass
                else:
                    if last_count != 0:
                        f.write(f'{chrom}\t{region_start+1}\t{index}\t{last_count}\n')
                    last_count = count
                    region_start = index
            else:
                if (count == last_count) and count != 0:
                    f.write(f'{chrom}\t{region_start+1}\t{index+1}\t{last_count}\n')

if __name__ == '__main__':
    main()
