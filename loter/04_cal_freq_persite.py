# -*- coding: utf-8 -*-
"""
Created on Fri Apr 30 22:50:24 2021

@author: Yudongcai

@Email: yudong_cai@163.com
"""

import gzip
import typer
import numpy as np
import pandas as pd
from collections import defaultdict



def main(ancfile: str = typer.Option(..., help="第一步跑出来的.anc.tsv.gz"),
         bagfile: str = typer.Option(..., help="第一步跑出来的.bag.tsv.gz"),
         outfile: str = typer.Option(..., help="输出文件.freq.tsv.gz"),
         popfile: str = typer.Option(..., help="样本分组信息,第一列个体ID(不是_1,_2的单倍体ID),第二列所属群体"),
         chunksize: int = typer.Option(10000, help="每次读取多少行"),
         cutoff: int = typer.Option(144, help="bag支持数的cutoff，默认160最大值，144就是90%"),
         refpops: str = typer.Option(..., help='参考群体ID，与第一步的这个参数一致')):
    
    # 样本信息
    group2haps = defaultdict(list)
    haps_used = []
    with open(popfile) as f:
        for line in f:
            tline = line.strip().split()
            haps_used.append(f'{tline[0]}_1')
            haps_used.append(f'{tline[0]}_2')
            group2haps[tline[1]].append(f'{tline[0]}_1')
            group2haps[tline[1]].append(f'{tline[0]}_2')
    group2nhaps = {}
    for group, haps in group2haps.items():
        group2nhaps[group] = len(haps)
    group_order = list(group2haps.keys())
    
    sources_order = refpops.split(',')
    source2index = {}
    for index, source in enumerate(sources_order):
        source2index[source] = index

    #输出header
    with gzip.open(outfile, 'wb') as f:
        outcols = ['chrom', 'pos']
        for source in sources_order:
            for group in group_order:
                name = f'{source}->{group}'
                outcols.append(name)
        outheader = '\t'.join(outcols) + '\n'
        f.write(outheader.encode())
        
    # 按需指定数据类型
    cols = ['chrom', 'pos', 'region'] + haps_used
    dtypes_anc = {'chrom': 'category', 'pos': np.uint32, 'region': 'category'}
    for hap in haps_used:
        dtypes_anc[hap] = np.int8
    dtypes_bag = {'chrom': 'category', 'pos': np.uint32, 'region': 'category'}
    for hap in haps_used:
        dtypes_bag[hap] = np.uint16
    
    # 分块处理
    with pd.read_csv(ancfile, usecols=cols, compression='gzip', sep='\t', dtype=dtypes_anc, chunksize=chunksize) as reader1, \
         pd.read_csv(bagfile, usecols=cols, compression='gzip', sep='\t', dtype=dtypes_bag, chunksize=chunksize) as reader2:
        for adf, bdf in zip(reader1, reader2):
            odf = adf[['chrom', 'pos']].copy()
            # 不满足cutoff的mask为-1
            adf[haps_used].values[bdf[haps_used].values < cutoff] = -1 
            for source in sources_order:
                source_index = source2index[source]
                for group in group_order:
                    name = f'{source}->{group}'
                    haps_group = group2haps[group]
                    nhaps = group2nhaps[group]
                    # 分group计算频率
                    odf[name] = np.sum(adf[haps_group].values == source_index, axis=1) / nhaps
            odf[outcols].to_csv(outfile, sep='\t', float_format='%.3f', mode='a', compression='gzip', header=None, index=False)


if __name__ == '__main__':
    typer.run(main)