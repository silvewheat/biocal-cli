# -*- coding: utf-8 -*-
"""
Created on Thu Apr 29 17:23:01 2021

@author: Yudongcai

@Email: yudong_cai@163.com
"""

import gzip
import pysam
import typer
import numpy as np
import pandas as pd


def open_bag(bagfile):
    with gzip.open(bagfile, 'rb') as f:
        tbxheader = np.array(f.readline().decode().split('\t'))
    tbx_bag = pysam.TabixFile(bagfile, parser=pysam.asTuple())
    return tbxheader, tbx_bag

def fetch_bag(tbx_bag, tbxheader: list, chrom: str, start: int, end: int):
    tdf = []
    for row in tbx_bag.fetch(chrom, start, end):
        tdf.append([int(x) for x in row[3:]])
    tdf = pd.DataFrame(tdf, columns=tbxheader[3:], dtype=np.uint8)
    return tdf

def load_indseg(indseg, haps: None):
    cols_used = ['chrom', 'start', 'end', 'nSNPs', 'hapID', 'source', 'group']
    dtypes = {'chrom': 'category', 'hapID': 'category',
              'source': 'category', 'group': 'category'}
    rdf = pd.read_csv(indseg, sep='\t', dtype=dtypes, usecols=cols_used)
    if haps:
        haps = {x.strip() for x in open(haps)}
        rdf = rdf.loc[rdf['hapID'].isin(haps), :]
    return rdf


def main(indseg: str = typer.Argument(..., help="合并后的每个hap的区间文件.indseg.tsv.gz"),
         bagfile: str = typer.Argument(..., help="bgzip压缩，tabix建索引的bag.tsv.bgz"),
         outfile: str = typer.Argument(..., help="输出文件名.indseg.addBag.tsv.gz"),
         haps: str = typer.Option(None, help='只输出这些单倍体，一行一个单倍体ID'),
         bagmax: int = typer.Option(160, help='bagging的最大值，不懂的话看文章')):
    """
    第一步的out.bag.tsv.gz需要进一步用bgzip压缩并建索引
    zcat out.bag.tsv.gz | bgzip -c > out.bag.tsv.bgz
    tabix -S 1 -s 1 -b 2 -e 2 -f out.bag.tsv.bgz
    """
    rdf = load_indseg(indseg, haps)
    tbxheader, tbx_bag = open_bag(bagfile)
    with gzip.open(outfile, 'wb') as f:
        cols_out = ['chrom', 'start', 'end', 'nSNPs', 'hapID', 
                    'source', 'group', 'meanBag']
        header = '\t'.join(cols_out) + '\n'
        f.write(header.encode())
    for name, group in rdf.groupby(['chrom', 'start', 'end']):
        chrom, start, end = name
        tdf = fetch_bag(tbx_bag, tbxheader, str(chrom), start, end)
        meanbag = tdf[group['hapID'].values].mean() / bagmax
        group['meanBag'] = group['hapID'].map(meanbag)
        group[cols_out].to_csv(outfile, sep='\t', header=False, 
                               float_format='%.3f', mode='a', index=False,
                               compression='gzip')
        

if __name__ == '__main__':
    typer.run(main)