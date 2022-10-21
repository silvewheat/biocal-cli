# -*- coding: utf-8 -*-
"""
Created on Tue Oct 16 21:40:08 2018

@author: YudongCai

@Email: yudongcai216@gmail.com
"""


import click
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA


def load_data(infile):
    """
    读纯数字的输入矩阵
    每行一个个体，每列一个特征
    """
    df = pd.read_csv(infile, sep='\t', header=None)
    return df.values


def pca(X):
    pca = PCA(n_components=10)
    pca.fit(X)
    components_ay = pca.transform(X)
    explained_variance_ratio = pca.explained_variance_ratio_
    return components_ay, explained_variance_ratio




@click.command()
@click.argument('infile')
@click.argument('outprefix')
def main(infile, outprefix):
    """
    输入的INFILE是纯数字矩阵，tab分割，每行一个个体，每列一个特征
    """
    X = load_data(infile)
    pcs, ratios = pca(X)
    pd.DataFrame(pcs).to_csv(f'{outprefix}.pcs', sep='\t', header=None, index=False)
    ratios.tofile(f'{outprefix}.explained_variance_ratio', sep='\n')



if __name__ == '__main__':
    main()