# -*- coding: utf-8 -*-
"""
Created on Sun Mar  8 10:38:45 2020

@author: YudongCai

@Email: yudongcai216@gmail.com
"""

import click
import numpy as np
import pandas as pd
import allel


def getAC(genotypes, variant_selection, sample_selection):
    return genotypes.subset(variant_selection, sample_selection).count_alleles()


@click.command()
@click.option('--vcffile')
@click.option('--pop1', help='sample ID list')
@click.option('--pop2', help='sample ID list')
@click.option('--binwidth', type=int, default=50000, help='滑动窗口大小')
@click.option('--stepsize', type=int, default=10000, help='滑动窗口步长')
@click.option('--outprefix')
def main(vcffile, pop1, pop2, binwidth, stepsize, outprefix):
    """
    计算pop1和pop2之间的Fst
    using the method of Hudson (1992) elaborated by Bhatia et al. (2013).
    """
    pop1 = [x.strip() for x in open(pop1)]
    pop2 = [x.strip() for x in open(pop2)]
    callset = allel.read_vcf(vcffile)
    allsamples = callset['samples']
    genotypes = allel.GenotypeChunkedArray(callset['calldata/GT'])
    variant_selection = np.full((genotypes.shape[0]+1), True) # 选择vcf中的全部位点
    sample_selection =  [True if x in pop1 else False for x in allsamples]
    ac1 = getAC(genotypes, variant_selection, sample_selection)
    sample_selection =  [True if x in pop2 else False for x in allsamples]
    ac2 = getAC(genotypes, variant_selection, sample_selection)
    num, den = allel.hudson_fst(ac1, ac2)
    fst = num / den
    meanFst = np.sum(num) / np.sum(den)
    print('meanFst: %s' % meanFst)
    chrom = callset['variants/CHROM']
    pos = callset['variants/POS']
    df = pd.DataFrame({'chrom': chrom, 'pos': pos, 'hudson_Fst': fst})
    df.to_csv(f'{outprefix}_persite.tsv.gz', sep='\t', index=False, na_rep='nan', compression='gzip')
    df['num'] = num
    df['den'] = den
    # sliding bins
    bdf = []
    for offset in range(0, binwidth, stepsize):
        df['bin_index'] = ((df['pos'].values - 1) - offset) // binwidth
        for group_name, gdf in df.groupby(by=['chrom', 'bin_index']):
            chrom, bin_index = group_name
            start = bin_index * binwidth + offset + 1
            if start < 0: # 开头几个窗口长度不足的就直接跳过
                continue
            end = start + binwidth - 1
            n_snp = gdf.shape[0]
            sum_num = gdf['num'].sum()
            sum_den = gdf['den'].sum()
            if sum_den > 0:
                meanFst = sum_num / sum_den
            else:
                meanFst = np.nan
            bdf.append([chrom, start, end, n_snp, meanFst])
    bdf = pd.DataFrame(bdf, columns=['chrom', 'start', 'end', 'n_snp', 'meanFst']).sort_values(by=['chrom', 'start'])
    bdf.to_csv(f'{outprefix}_meanFst.tsv.gz', index=False, compression='gzip', sep='\t', float_format='%.3f')



if __name__ == '__main__':
    main()


