import typer
import allel
import numpy as np
import pandas as pd


def main(vcffile: str = typer.Argument(..., help="输入的vcf文件"),
         outfile: str = typer.Argument(..., help="输出文件名(tsv文件)"),
         samples: str = typer.Option(None, help="仅使用这些样本（一行一个个体）")):
    """
    统计指定VCF文件中每个个体的测序深度分布，vcf FMT中需要包含DP字段。
    """
    if samples:
        samples = [x.strip() for x in open(samples)]
    callset = allel.read_vcf(vcffile, samples=samples, fields=['samples', 'calldata/DP'])
    samples = callset['samples']
    df = pd.DataFrame(callset['calldata/DP'].sum(axis=2), columns=callset['samples'])
    df = df.replace({-1: 0})
    nsites, nsamples = df.shape
    print(f'{nsamples} samples, {nsites} sites loaded.')
    cdf = {}
    cdf['sample'] = df.columns
    cdf['1X_all'] = np.sum(df.values >= 1, axis=0) / nsites
    cdf['3X_all'] = np.sum(df.values >= 3, axis=0) / nsites
    cdf['5X_all'] = np.sum(df.values >= 5, axis=0) / nsites
    cdf['10X_all'] = np.sum(df.values >= 10, axis=0) / nsites
    cdf['20X_all'] = np.sum(df.values >= 20, axis=0) / nsites
    cdf = pd.DataFrame(cdf)
    cdf.to_csv(outfile, sep='\t', index=False)


if __name__ == '__main__':
    typer.run(main)