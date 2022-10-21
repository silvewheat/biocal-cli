# -*- coding: utf-8 -*-
"""
Created on Sat May  1 14:38:31 2021

@author: Yudongcai

@Email: yudong_cai@163.com
"""

import typer
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt



def plot(bdf, chrom, start, end, binwidth, source2color, group_order, sources_order, outfile):
    rdf = bdf.loc[(bdf['binleft']>=start)&(bdf['binleft']<=(end-binwidth)), :]
    fig, ax = plt.subplots(len(group_order), 1, figsize=(20, 5), sharey=True, sharex=True)
    for ngroup, group in enumerate(group_order):
        botoom = np.zeros(rdf.shape[0])
        for nsource, source in enumerate(sources_order):
            ax[ngroup].set_ylabel(group)
            ax[ngroup].bar(x=rdf['binleft'].values,
                   height=rdf[f'{source}->{group}'].values,
                   width=binwidth,
                   bottom=botoom,
                   color=source2color[source],
                   linewidth=0,
                   align='edge',
                   label=source)
            ax[ngroup].hlines(y=0.5, ls='dashed', xmin=start, xmax=end, color='red', lw=0.5)
            botoom += rdf[f'{source}->{group}'].values
    ax[-1].set_ylim([0, 1])
    ax[-1].set_xlim([rdf['binleft'].min(), rdf['binleft'].max() + binwidth])
    plt.legend(bbox_to_anchor=(1,0), loc="lower left", ncol=1)
    ax[0].set_title(f'chr{chrom}:{start//1000000}M-{end//1000000}M')
    plt.savefig(outfile, dpi=200)
    plt.close()


def load_freq(freqfile, group_order, sources_order):
    cols = ['chrom', 'pos']
    dtypes = {'chrom': 'category', 'pos': np.uint32}
    for source in sources_order:
        for group in group_order:
            comb = f'{source}->{group}'
            dtypes[comb] = np.float16
            cols.append(comb)
    df = pd.read_csv(freqfile, sep='\t', dtype=dtypes)
    return df


def process(cdf, binwidth):
    cdf['binindex'] = cdf['pos'].values // binwidth
    bdf = cdf.groupby('binindex').mean().reset_index()
    bdf['binleft'] = bdf['binindex'] * binwidth
    bdf['binleft'] = bdf['binleft'].astype(int)
    return bdf


def main(freqfile: str = typer.Option(..., help="第四步跑出来的.freq.tsv.gz"),
         grouporder: str = typer.Option(..., help="画图时query群体的顺序，逗号分割"),
         sourcesorder: str = typer.Option(..., help="画图时source群体的顺序，逗号分割"),
         sourcecolor: str = typer.Option(..., help="source的颜色，两列，一列ID一列十六进制颜色(#开头)"),
         binwidth: int = typer.Option(..., help="窗口大小"),
         outprefix: str = typer.Option(..., help="输出文件前缀")):
    
    group_order = grouporder.strip().split(',')
    sources_order = sourcesorder.strip().split(',')
    
    source2color = {x.split()[0]: x.strip().split()[1] for x in open(sourcecolor)}
    
    df = load_freq(freqfile, group_order, sources_order)
    
    for chrom in df['chrom'].unique():
        cdf = df.loc[df['chrom']==chrom, :].copy()
        bdf = process(cdf, binwidth)
        for start in range(1, bdf['binleft'].max() + binwidth, 10_000_000): # 每个图片画10M
            end = start + 9_999_999
            print(chrom, start, end)
            outfile = f'{outprefix}_chr{chrom}_{end//1000000}M.jpg'
            plot(bdf, chrom, start, end, binwidth, source2color, group_order, sources_order, outfile)


if __name__ == '__main__':
    typer.run(main)