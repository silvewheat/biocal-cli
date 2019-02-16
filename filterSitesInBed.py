# -*- coding: utf-8 -*-
"""
Created on Tue Oct 16 15:12:22 2018

@author: YudongCai

@Email: yudongcai216@gmail.com
"""

import os
import gzip
import click
import numpy as np
import pandas as pd



def count_header(gzvcf):
    with gzip.open(gzvcf, 'rb') as f:
        for nline, line in enumerate(f):
            if line.decode()[0] != '#':
                return nline


@click.command()
@click.option('--sitesfile', help='Chrom\\tPos', default=None)
@click.option('--gzvcf', help='bgzip vcffile', default=None)
@click.option('--bedfile', help='Chrom\\tStart\\tEnd')
@click.option('--outsites', help='output file')
def main(sitesfile, gzvcf, bedfile, outsites):
    """
    output sites (in sitesfile or gzvcffile) within the bedfile
    """
    try:
        os.remove(outsites)
    except Exception:
        pass
    if gzvcf:
        nheader = count_header(gzvcf)
        sdf = pd.read_csv(gzvcf, skiprows=nheader, header=None, names=['chrom','pos'],
                         usecols=[0,1], sep='\t', dtype={'chrom': str, 'pos': int})
    elif sitesfile:
        sdf = pd.read_csv(sitesfile, sep='\t', header=None, names=['chrom', 'pos'],
                          dtype={'chrom': str, 'pos': int})
    bdf = pd.read_csv(bedfile, sep='\t', header=None, names=['chrom', 'start', 'end'],
                      dtype={'chrom': str, 'start': int, 'end': int})
    for chrom, start, end in bdf[['chrom', 'start', 'end']].values:
        sdf.loc[(sdf['chrom']==chrom) & (sdf['pos']>=start) & (sdf['pos']<=end), :]\
            .to_csv(outsites, mode='a', header=False, sep='\t', index=False)


if __name__ == '__main__':
    main()
