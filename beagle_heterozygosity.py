# -*- coding: utf-8 -*-
"""
Created on Fri May 10 08:58:27 2019

@author: YudongCai

@Email: yudongcai216@gmail.com
"""

import gzip
import click
import numpy as np
import pandas as pd





@click.command()
@click.option('--beaglefile', help='压缩后的beagle文件')
@click.option('--chuncksize', help='一次处理XX行，默认10000', default=10000, type=int)
@click.option('--outfile', help='输出压缩后文件的名字(不压缩)')
def main(beaglefile, chuncksize, outfile):
    """
    从bealge.gz文件中计算杂合度Individual heterozygosity
    方法参考[1] A. Fages et al., “Tracking Five Millennia of Horse Management with Extensive Ancient Genome Time Series.,” Cell, vol. 0, no. 0, Apr. 2019.
    beagle格式：http://www.popgen.dk/angsd/index.php/Genotype_Likelihoods#Beagle_format
    allele codes as 0=A, 1=C, 2=G, 3=T
    Genotype likelihood for the major/major, major/minor, minor/minor
    """
    with gzip.open(beaglefile, 'rb') as f1, gzip.open(outfile, 'w') as f2:
        header = f1.readline().decode().strip().split()
        samples = header[3::3]
        outheader = 'SMID\tHet\tHet_rmTransi\n'
        f2.write(outheader.encode())
    ucols = [0, 1, 2] + list(range(4, len(header), 3))
    names = ['marker', 'allele1', 'allele2']
    reader = pd.read_csv(beaglefile, sep='\t', skiprows=1, header=None, names=names,
                         na_values='0.333333', iterator=True, usecols=ucols)
    loop = True
    while loop:
        try:
            chunk = reader.get_chunk(chuncksize)
            df1 = chunk[samples].sum(axis=0) / chunk[samples].count(axis=0)
            df2 = chunk[samples].count(axis=0)
            odf1 = pd.concat([df1, df2], axis=1, keys=['Het_all', 'Count_all'])

            transversion = np.fabs(chunk['allele1'] - chunk['allele2'])!=2
            df3 = chunk.loc[transversion, samples].sum(axis=0) / chunk.loc[transversion, samples].count(axis=0)
            df4 = chunk.loc[transversion, samples].count(axis=0)
            odf2 = pd.concat([df3, df4], axis=1, keys=['Het_transversion', 'Count_transversion'])

            odf =

            if rmtansi:
                chunk = chunk.loc[np.fabs(chunk[1] - chunk[2])!=2, :]
            chunk.to_csv(outfile, sep='\t', index=False, header=False, compression='gzip', mode='ab')
        except StopIteration:
            loop = False
            print('Done!')

