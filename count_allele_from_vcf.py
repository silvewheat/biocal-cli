# -*- coding: utf-8 -*-
"""
Created on Tue Aug 10 17:35:37 2021

@author: Yudongcai

@Email: yudong_cai@163.com
"""
import click
import allel
import pandas as pd
from collections import defaultdict




@click.command()
@click.option('--vcffile')
@click.option('--groupfile', help='sampleID\tgroupID')
@click.option('--outfile', help='outfile name')
def main(vcffile, groupfile, outfile):
    group2inds = defaultdict(list)
    with open(groupfile) as f:
        for line in f:
            sampleID, groupID = line.strip().split()
            group2inds[groupID].append(sampleID)
    callset = allel.read_vcf(vcffile, fields=['variants/CHROM', 'variants/POS', 'variants/REF', 'variants/ALT'],
                             numbers={'ALT': 1}) # 第2个及以上的ALT将被忽略（但是位点还在） 多等位推荐把不同ALT分开在vcf的不同行
    df = pd.DataFrame({'chr': callset['variants/CHROM'],
                       'pos': callset['variants/POS'],
                       'REF': callset['variants/REF'],
                       'ALT': callset['variants/ALT']})
    for group, samples in group2inds.items():
        print(group)
        print(samples)
        callset = allel.read_vcf(vcffile, samples=samples,
                                 fields=['samples', 'calldata/GT'])
        GT_ay = allel.GenotypeArray(callset['calldata/GT'])
        af = GT_ay.count_alleles().to_frequencies()
        df[f'{group}_reffreq'] = af[:, 0]
        if af.shape[1] > 1: # ALT如果频率都是0的话，就只会有一列REF的频率了
            df[f'{group}_altfreq'] = af[:, 1] # 第一个ALT的频率
        else:
            df[f'{group}_altfreq'] = .0
        # 基因型计数
        df[f'{group}_nCalled'] = GT_ay.count_called(axis=1)
        df[f'{group}_nHomRef'] = GT_ay.count_hom_ref(axis=1)
        df[f'{group}_nHomAlt'] = GT_ay.count_hom_alt(axis=1)
        df[f'{group}_nHet'] = GT_ay.count_het(axis=1)


    df.to_csv(outfile, sep='\t', index=False, float_format='%.3f', na_rep='nan')


if __name__ == '__main__':
    main()

