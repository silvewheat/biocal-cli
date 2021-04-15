# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 21:47:54 2021

@author: Yudongcai

@Email: yudong_cai@163.com
"""

import gzip
import typer
import allel
import numpy as np
import pandas as pd
from collections import defaultdict
import loter.locanc.local_ancestry as lc



def vcf2npy(vcffile, samples, region=None, return_more=False):
    fields = ['samples', 'calldata/GT']
    if return_more:
        fields=['samples', 'calldata/GT', 'variants/CHROM', 'variants/POS']
    callset = allel.read_vcf(vcffile, samples=samples, region=region,
                             fields=fields)
    haplotypes_1 = callset['calldata/GT'][:,:,0]
    haplotypes_2 = callset['calldata/GT'][:,:,1]

    m, n = haplotypes_1.shape
    mat_haplo = np.empty((2*n, m))
    mat_haplo[::2] = haplotypes_1.T
    mat_haplo[1::2] = haplotypes_2.T

    keep_samples = callset['samples']
    
    if return_more:
        return mat_haplo.astype(np.uint8), keep_samples, callset['variants/CHROM'], callset['variants/POS']
    else:
        return mat_haplo.astype(np.uint8)


def load_group(groupfile):
    sample2group = {}
    group2samples = defaultdict(list)
    with open(groupfile) as f:
        for line in f:
            sample, group = line.strip().split()
            sample2group[sample] = group
            group2samples[group].append(sample)
    return sample2group, group2samples


def load_regions(regionfile):
    regions = []
    with open(regionfile) as f:
        for line in f:
            tline = line.strip().split()
            regions.append(f'{tline[0]}:{tline[1]}-{tline[2]}')
    return regions


def main(vcffile: str = typer.Argument(..., help="总vcf文件"),
         outfile: str = typer.Argument(..., help='输出文件名(gzipped), XX.tsv.gz'),
         groupfile: str = typer.Option(..., help='样本分群信息，两列，一列vcf中的样本ID，一列对应的群体ID'),
         regionfile: str = typer.Option(..., help='区域文件，三列，chrom\\tstart\\tend'),
         refpops: str = typer.Option(..., help='参考群体ID，此ID须存在于groupfile的第二列中。多个群体用,分割'),
         querypops: str = typer.Option(..., help='待检测群体ID，此ID须存在于groupfile的第二列中。多个群体用,分割'),
         threads: int = typer.Option(1, help='使用的线程数'),
         nbags: int = typer.Option(20, help='number of resampling in the bagging')):
    
    sample2group, group2samples = load_group(groupfile)
    
    querysamples = []
    for group in querypops.split(','):
        querysamples.extend(group2samples[group])
        
    hapIDs = [[f'{i}_1', f'{i}_2'] for i in querysamples]
    hapIDs = [i for j in hapIDs for i in j]
    hapID2group = {x: sample2group[x[:-2]] for x in hapIDs}
    
    cols = ['chrom', 'pos', 'region', 'hapID', 'group', 'sourceIndex', 'nb_bagging']
    with gzip.open(outfile, 'wb') as f:
        header = '\t'.join(cols) + '\n'
        f.write(header.encode())

    # 分区域计算
    regions = load_regions(regionfile)
    for region in regions:
        print(region)

        H_query, query_samples, chrs, sites = vcf2npy(vcffile, querysamples, region=region, return_more=True)

        H_refs = []
        for group in refpops.split(','):
            H_ref = vcf2npy(vcffile, group2samples[group], region=region)
            H_refs.append(H_ref)

        print(f'{len(sites)} sites loaded.')
        
        res_loter = lc.loter_local_ancestry(H_refs, H_query, num_threads=threads, nb_bagging=nbags)
        
        df = pd.DataFrame(res_loter[0].T, columns=hapIDs)
        df['chrom'] = chrs
        df['pos'] = sites
        mdf = df.melt(id_vars=['chrom', 'pos'], value_vars=hapIDs, var_name='hapID', value_name='sourceIndex')
        mdf['group'] = mdf['hapID'].map(hapID2group)

        bdf = pd.DataFrame(res_loter[1].T, columns=hapIDs)
        bdf['chrom'] = chrs
        bdf['pos'] = sites
        mbdf = bdf.melt(id_vars=['chrom', 'pos'], value_vars=hapIDs, var_name='hapID', value_name='nb_bagging')
        mbdf['group'] = mbdf['hapID'].map(hapID2group)
        
        mdf = pd.merge(mdf, mbdf, on=['chrom', 'pos', 'hapID', 'group'], how='inner')
        mdf['region'] = region # 这个很关键，因为是分区域算的，区域之间不连续，后续位点合并为区间只能在一个区域内进行
        mdf[cols].to_csv(outfile, sep='\t', index=False, header=None, mode='a', compression='gzip')

if __name__ == '__main__':
    typer.run(main)

# 添加区间去冗余的参数
# 去掉无变异位点的参数