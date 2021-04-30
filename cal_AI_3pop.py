# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 15:35:14 2020

@author: YudongCai

@Email: yudongcai216@gmail.com
"""

import allel
import click
import numpy as np
import pandas as pd


def select_samples(samples_all, samples_query):
    """
    从samples_all里面筛选samples_query，产生布尔掩码
    """
    return [True if sample in samples_query else False for sample in samples_all]

def load_vcf2array(vcffile, samples_queried, outsamples: None):
    callset = allel.read_vcf(vcffile, samples=samples_queried,
                             fields=['samples', 'calldata/GT',  'variants/CHROM', 'variants/POS'])
    gt_array = callset['calldata/GT'] # 三维array
    samples_all = callset['samples']
    pos_array = callset['variants/POS']
    chrom_array = callset['variants/CHROM']
    
    # 样本是否都在vcf里
    if len(samples_queried) != len(samples_all):
        samples_notexist = set(samples_queried) - set(samples_all)
        print(f'{len(samples_notexist)} samples not exist in the vcf file:')
        print(', '.join(samples_notexist))
    
    # 只保留双等位, 多等位以后再考虑
    selection_biallelic = np.max(np.max(gt_array, axis=2), axis=1) < 2
    gt_array = gt_array[selection_biallelic, :, :]
    pos_array = pos_array[selection_biallelic]
    chrom_array = chrom_array[selection_biallelic]
    n_sites, n_samples, n_hap = gt_array.shape
    print(f'{n_sites} biallelic sites were remained.')
    
    # 把outsamples中最高频率的allele设置为alt
    if outsamples:
        selection_outsamples = select_samples(samples_all, outsamples)
        print(f'{np.sum(selection_outsamples)} outgroup samples in vcf file.')
        gt_array_out = gt_array[:, selection_outsamples, :].reshape(n_sites, np.sum(selection_outsamples)*n_hap)
        selection_swtich = np.sum(gt_array_out==0, axis=1) > np.sum(gt_array_out>0, axis=1) # ref(0)的数量比非ref(!=0)但不是miss(-1)的数量多
        print(f'Swtich REF and ALT in {np.sum(selection_swtich)} sites.')
        assert gt_array.min() >= -1
        assert gt_array.max() <= 1
        gt_swtich = gt_array[selection_swtich, :, :]
        gt_swtich[gt_swtich == 1] = 9
        gt_swtich[gt_swtich == 0] = 1
        gt_swtich[gt_swtich == 9] = 0
        gt_array[selection_swtich, :, :] = gt_swtich
    return gt_array, callset['samples'], pos_array, chrom_array

def cal_alt1_freq(gt_array):
    """
    gt_array为load_vcf2array产生的3维ndarray
    注意，返回的是第一个ALT的frequency
    """
    return allel.GenotypeArray(gt_array).count_alleles().to_frequencies()[:, 1]


@click.command()
@click.option('--vcffile', help='输入的vcf文件')
@click.option('--popa', help='群体A(非渗入群体)的ID列表，一行一个')
@click.option('--popb', help='群体B(被渗入群体)的ID列表，一行一个')
@click.option('--popc', help='群体C(渗入来源群体)的ID列表，一行一个')
@click.option('--binwidth', type=int, default=50000, help='滑动窗口大小')
@click.option('--stepsize', type=int, default=10000, help='滑动窗口步长')
@click.option('--outprefix', help='输出文件前缀')
def main(vcffile, popa, popb, popc, binwidth, stepsize, outprefix):
    """
    U_A,B,C(w,x,y)
    A是非渗入群体，B是被渗入群体，C是渗入来源群体
    在窗口内A中频率小于w，B中大于x，C中大于y的SNP位点数即为U_A,B,C(w,x,y)
    详见：Signatures of Archaic Adaptive Introgression in Present-Day Human Populations
    """
    samples_popA = [x.strip() for x in open(popa)]
    samples_popB = [x.strip() for x in open(popb)]
    samples_popC = [x.strip() for x in open(popc)]
    samples_queried = samples_popA + samples_popB + samples_popC
    assert len(set(samples_queried)) == len(samples_queried), "样本ID有重复"
    gt_array, samples_all, pos_array, chrom_array = load_vcf2array(vcffile, samples_queried, outsamples=samples_popC)
    
    # 群体C频率
    selection_popC = select_samples(samples_all, samples_popC)
    af_popC = cal_alt1_freq(gt_array[:, selection_popC, :]) # 这个返回的是第一个ALT的frequency
    # 只保留群体C中频率百分百的位点, ALT按popC转换过了，所以直接算ALT的频率
    selection_freq_popC = af_popC == 1
    gt_array = gt_array[selection_freq_popC, :, :]
    pos_array = pos_array[selection_freq_popC]
    chrom_array = chrom_array[selection_freq_popC]
    af_popC = af_popC[selection_freq_popC]

    # 群体AB的频率，使用筛选后的位点计算
    selection_popA = select_samples(samples_all, samples_popA)
    af_popA = cal_alt1_freq(gt_array[:, selection_popA, :]) # 这个返回的是第一个ALT的frequency
    selection_popB = select_samples(samples_all, samples_popB)
    af_popB = cal_alt1_freq(gt_array[:, selection_popB, :]) # 这个返回的是第一个ALT的frequency

    # 先把频率文件保存一下
    df = pd.DataFrame({'chrom': chrom_array, 'pos': pos_array, 'popA': af_popA, 'popB': af_popB, 'popC': af_popC})
    df.to_csv(f'{outprefix}_altFreq.tsv.gz', index=False, compression='gzip', sep='\t', float_format='%.3f')
    print('freq file saved.')
    
    # 过滤满足指标要求的位点
    selection_popA_1percent = df['popA'].values < 0.01
    df = df.iloc[selection_popA_1percent, :]
    
    selection_popC_1percent = df['popC'].values == 1
    df = df.iloc[selection_popC_1percent, :]
    print(f'{df.shape[0]} sites remained after frequency filtering.')
    
    # 滑动窗口计算统计量
    odf = []
    for offset in range(0, binwidth, stepsize):
        df['bin_index'] = ((df['pos'].values - 1) - offset) // binwidth
        for group_name, gdf in df.groupby(by=['chrom', 'bin_index']):
            chrom, bin_index = group_name
            start = bin_index * binwidth + offset + 1
            if start < 0: # 开头几个窗口长度不足的就直接跳过
                continue
            end = start + binwidth - 1
            n_snp = gdf.shape[0]
            print(chrom, start, end, n_snp)
#            Q_1_100_q90, Q_1_100_q95, Q_1_100_q100 = np.quantile(gdf['popB'].values, [0.9, 0.95, 1]) # A和B群体的频率已经提前过滤了 New in version 1.15.0.
            Q_1_100_q90, Q_1_100_q95, Q_1_100_q100 = np.percentile(gdf['popB'].values, [90, 95, 100])
            U_1_10_100 = np.sum(gdf['popB'].values >= 0.1)
            U_1_20_100 = np.sum(gdf['popB'].values >= 0.2)
            U_1_50_100 = np.sum(gdf['popB'].values >= 0.5)
            U_1_80_100 = np.sum(gdf['popB'].values >= 0.8)
            odf.append([chrom, start, end, n_snp, Q_1_100_q90, Q_1_100_q95, Q_1_100_q100, U_1_10_100, U_1_20_100, U_1_50_100, U_1_80_100])

    odf = pd.DataFrame(odf, columns=['chrom', 'start', 'end', 'n_snp', 'Q90', 'Q95', 'Q100', 'U10', 'U20', 'U50', 'U80']).sort_values(by=['chrom', 'start'])
    odf.to_csv(f'{outprefix}_stat.tsv.gz', index=False, compression='gzip', sep='\t', float_format='%.3f')

if __name__ == '__main__':
    main()
