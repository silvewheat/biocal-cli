# -*- coding: utf-8 -*-
"""
Created on Wed Aug 25 15:31:28 2021

@author: Yudongcai

@Email: yudong_cai@163.com
"""

import typer
import numpy as np
import pandas as pd


def load_paf(paffile):
    df = pd.read_csv(paffile, sep='\t', header=None, usecols=range(12),
                 names=['qname', 'qlen', 'qstart', 'qend', 'strand', 'tname', 'tlen', 'tstart', 'tend', 'residue', 'alnlen', 'mquality'])
    return df


def main(paffile: str = typer.Argument(..., help="minimap2输出的paf格式比对结果"),
         outfile: str = typer.Argument(..., help="输出的每条query的coverage统计结果"),
         minqcov: float = typer.Option(0.1, help="query覆盖的的最小值，低于这个不输出")):
    """
    从minimap2的比对结果中进行汇总，统计每条query序列在每个target序列中的覆盖率（qcov）。
    如果与某个target序列的覆盖率低于minqcov则不输出。
    """
    df = load_paf(paffile)
    sdf = df.groupby(['qname', 'tname']).agg({'alnlen': 'sum', 'qlen': 'max', 'tlen': 'max'}).reset_index()
    sdf['qcov'] = sdf['alnlen'].values / sdf['qlen'].values
    sdf['tcov'] = sdf['alnlen'].values / sdf['tlen'].values
    sdf = sdf.sort_values(['qlen', 'qcov'], ascending=False)
    sdf.loc[sdf['qcov']>=minqcov, :].to_csv(outfile, sep='\t', index=False)

if __name__ == '__main__':
    typer.run(main)
    
    
    
    