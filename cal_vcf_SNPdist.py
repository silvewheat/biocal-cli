import typer
import allel
import numpy as np
import pandas as pd
from itertools import combinations


def cal_dist(df, sample1, sample2):
    tdf = df[[sample1, sample2]].dropna()
    nsites = tdf.shape[0]
    if nsites > 0:
        dist = tdf.diff(axis=1)[sample2].abs().sum() / (nsites * 2) # 按二倍体去计算的 ！！！！多等位部分其实不能这么算！！！！！
    else:
        dist = np.nan
    return dist, nsites


def main(vcffile: str = typer.Argument(..., help="输入的vcf文件"),
         outfile: str = typer.Argument(..., help="输出文件名(tsv文件)"),
         samples: str = typer.Option(None, help="仅使用这些样本（一行一个个体）")):
    """
    在指定的vcf文件中计算指定的samples两两*二倍体*之间的差异SNP百分比（针对*二等位*）
    """
    if samples:
        samples = [x.strip() for x in open(samples)]
    callset = allel.read_vcf(vcffile, samples=samples, fields=['samples', 'calldata/GT'])
    samples = callset['samples']
    df = pd.DataFrame(callset['calldata/GT'].sum(axis=2), columns=callset['samples'])
    df = df.replace({-2: np.nan})
    print(f'{len(samples)} samples, {df.shape[0]} sites loaded.')
    odf = []
    for sample1, sample2 in combinations(samples, 2):
        print(sample1, sample2)
        dist, nsites = cal_dist(df, sample1, sample2)
        odf.append([sample1, sample2, nsites, dist])
        odf.append([sample2, sample1, nsites, dist])
    odf = pd.DataFrame(odf, columns=['sample1', 'sample2', 'nsites', 'dist']).sort_values(by=['sample1', 'sample2'])
    odf.to_csv(outfile, sep='\t', index=False)

if __name__ == '__main__':
    typer.run(main)
