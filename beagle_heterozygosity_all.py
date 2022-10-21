# -*- coding: utf-8 -*-
"""
Created on Fri May 10 09:50:13 2019

@author: YudongCai

@Email: yudongcai216@gmail.com
"""


import gzip
import click
import numpy as np
import pandas as pd





@click.command()
@click.option('--beaglefile', help='压缩后的beagle文件')
@click.option('--outfile', help='输出压缩后文件的名字(不压缩)')
def main(beaglefile, outfile):
    """
    从bealge.gz文件中计算杂合度Individual heterozygosity
    方法参考[1] A. Fages et al., “Tracking Five Millennia of Horse Management with Extensive Ancient Genome Time Series.,” Cell, vol. 0, no. 0, Apr. 2019.
    beagle格式：http://www.popgen.dk/angsd/index.php/Genotype_Likelihoods#Beagle_format
    allele codes as 0=A, 1=C, 2=G, 3=T
    Genotype likelihood for the major/major, major/minor, minor/minor
    """
    print(__doc__)
    with gzip.open(beaglefile, 'rb') as f1:
        header = f1.readline().decode().strip().split()
        samples = header[3::3]
    ucols = [0, 1, 2] + list(range(4, len(header), 3))
    names = ['marker', 'allele1', 'allele2'] + samples
    df = pd.read_csv(beaglefile, sep='\t', skiprows=1, header=None, names=names,
                         na_values='0.333333', usecols=ucols)
    df1 = df[samples].sum(axis=0) / df[samples].count(axis=0)
    df2 = df[samples].count(axis=0)
    odf1 = pd.concat([df1, df2], axis=1, keys=['Het_all', 'Count_all'])

    transversion = np.fabs(df['allele1'] - df['allele2'])!=2
    df3 = df.loc[transversion, samples].sum(axis=0) / df.loc[transversion, samples].count(axis=0)
    df4 = df.loc[transversion, samples].count(axis=0)
    odf2 = pd.concat([df3, df4], axis=1, keys=['Het_transversion', 'Count_transversion'])

    odf = pd.concat([odf1, odf2], axis=1).reset_index()
    odf.columns = ['SMID', 'Het_all', 'Count_all', 'Het_transversion', 'Count_transversion']
    odf.to_csv(outfile, index=False, sep='\t', float_format='%.6f')


if __name__ == '__main__':
    main()


