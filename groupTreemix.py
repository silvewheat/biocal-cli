# -*- coding: utf-8 -*-
"""
Created on Sun Feb 10 15:14:08 2019

@author: YudongCai

@Email: yudongcai216@gmail.com
"""

import click
import numpy as np
import pandas as pd
from collections import defaultdict
from collections import OrderedDict


@click.command()
@click.option('--treemix', help='treemix input file')
@click.option('--groupfile', help='two cols, sampleID groupID')
@click.option('--outfile', help='output file name')
def main(treemix, groupfile, outfile):
    """
    merge treemix input file according to thier group (--groupfile)
    """
    groups = defaultdict(list)
    with open(groupfile) as f:
        for line in f:
            sample, group = line.strip().split()
            groups[group].append(sample)
    df = pd.read_csv(treemix, sep='\s+', header=0)
    first_haplotype = df.iloc[:, :].applymap(lambda s: int(s.split(',')[0]))
    second_haplotype = df.iloc[:, :].applymap(lambda s: int(s.split(',')[1]))
    fdf = OrderedDict()
    sdf = OrderedDict()
    for group, samples in groups.items():
        fdf[group] = first_haplotype[samples].sum(axis=1)
        sdf[group] = second_haplotype[samples].sum(axis=1)
    fdf = pd.DataFrame(fdf)
    sdf = pd.DataFrame(sdf)
    mdf = fdf.astype('str') + ',' + sdf.astype('str')
    if outfile.split('.')[-1] == 'gz':
        comp = 'gzip'
    else:
        comp = None
    mdf.to_csv(outfile, sep=' ', index=False, compression=comp)


if __name__ == '__main__':
    main()


