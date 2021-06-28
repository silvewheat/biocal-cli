# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 15:21:02 2021

@author: Yudongcai

@Email: yudong_cai@163.com
"""

import typer
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns


def load_dp(dpfile):
    df = pd.read_csv(dpfile, sep='\t', header=None,
                     names=['contig', 'start', 'end', 'depth']).dropna()
    return df, df['depth'].values

def esti_peak(dp, low_cutoff):
    """
    估算正常测序深度的峰值，只保留大于low_cutoff的数据，为避免大量低深度序列带来的错误（低深度的计数超过正常深度的计数）
    只保留一位小数，这对于取众数很重要！！！
    """
    modes, counts = stats.mode([round(x, 1) for x in dp[dp>5]])
    return modes[0], counts[0]

def esti_std(dp, peak):
    """
    使用大于等于peak的数据对称后得到的正态分布来计算std
    """
    dp2 = []
    dp2.extend(dp[dp>=peak])
    dp2.extend((dp[dp>peak] - dp.max())*-1 - dp.max() + peak*2)
    std = np.std(dp2)
    # ± 2SD
    upper_boundary = peak + std * 2
    lower_boundary = peak - std * 2
    return std, upper_boundary, lower_boundary, dp2

def plot_raw_dp(dp, peak, peak_count, std, upper_boundary, lower_boundary, outprefix):
    fig, ax = plt.subplots(1,1, figsize=(10,5))
    sns.histplot(dp, ax=ax, binwidth=0.1)
    y_max = peak_count * 1.1
    x_max = peak + std * 3
    ax.vlines(peak, 0, y_max, color='r', ls='dashed')
    ax.vlines(upper_boundary, 0, y_max, color='b', ls='dashed')
    ax.vlines(lower_boundary, 0, y_max, color='g', ls='dashed')
    ax.set_xlim([0, x_max])
    ax.set_ylim([0, y_max])
    plt.savefig(f'{outprefix}_rawDP.jpg', dpi=200)
    plt.close()

def plot_dp2(dp2, peak, peak_count, std, upper_boundary, lower_boundary, outprefix):
    fig, ax = plt.subplots(1,1, figsize=(10,5))
    sns.histplot(dp2, ax=ax, binwidth=0.1)
    y_max = peak_count * 1.1
    x_max = peak + std * 3
    ax.vlines(peak, 0, y_max, color='r', ls='dashed')
    ax.vlines(upper_boundary, 0, y_max, color='b', ls='dashed')
    ax.vlines(lower_boundary, 0, y_max, color='g', ls='dashed')
    ax.set_xlim([0, x_max])
    ax.set_ylim([0, y_max])
    plt.savefig(f'{outprefix}_DP2.jpg', dpi=200)
    plt.close()

def save_results(df, peak, std, upper_boundary, lower_boundary, outprefix):
    # 保存cutoff
    with open(f'{outprefix}.cutoff.tsv', 'w') as f:
        f.write('peak\tstd\tupper_boundary\tlower_boundary\n')
        f.write(f'{peak:.1f}\t{std:.2f}\t{upper_boundary:.2f}\t{lower_boundary:.2f}\n')
    # 保存满足标准的窗口ID
    tdf = df.loc[(df['depth']<upper_boundary)&(df['depth']>=lower_boundary), ['contig', 'start', 'end']]
    (tdf['contig'] + ':' + tdf['start'].astype(str) + '-' + tdf['end'].astype(str)).to_csv(f'{outprefix}.pass.tsv', index=False, header=False)

def main(dpfile: str = typer.Option(..., help="mosdepth生成的滑动窗口.regions.bed.gz"),
         outprefix: str = typer.Option(..., help="输出文件前缀"),
         low_cutoff: float = typer.Option(5, help="估算peak时只保留大于low_cutoff的数据")):
    """
    mosdepth统计BAM文件深度，使用--by 1000生成1K滑动窗口的深度统计结果.regions.bed.gz(也可使用其它窗口大小)。
    之后使用该脚本估算基因组正常区域的深度分布峰值，以及正常区域的深度分布范围。
    """
    df, dp = load_dp(dpfile)
    peak, peak_count = esti_peak(dp, low_cutoff)
    std, upper_boundary, lower_boundary, dp2 = esti_std(dp, peak)
    plot_raw_dp(dp, peak, peak_count, std, upper_boundary, lower_boundary, outprefix)
    plot_dp2(dp2, peak, peak_count, std, upper_boundary, lower_boundary, outprefix)
    save_results(df, peak, std, upper_boundary, lower_boundary, outprefix)


if __name__ == '__main__':
    typer.run(main)