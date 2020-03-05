# -*- coding: utf-8 -*-
"""
Created on Thu Jan  2 09:49:45 2020

@author: YudongCai

@Email: yudongcai216@gmail.com
"""

import click
import numpy as np
import pandas as pd
from itertools import cycle



def getJackKnife(df, ABBAname,BABAname,BBAAname):
    ## 过滤没有点的block
    print(df.shape)
    print('过滤6:')
    colWeights = df.iloc[:, 6:].sum(axis=1).values
    zeroIdx = colWeights != 0
    df = df.loc[zeroIdx, :]
    colWeights = colWeights[zeroIdx]
    print(df.shape)

    # 过滤分母是0的block
    print('过滤 colDen')
    colDen = df[ABBAname].sum(axis=1).values + df[BABAname].sum(axis=1).values
    zeroIdx = colDen!=0
    df = df.loc[zeroIdx, :]
    colWeights = colWeights[zeroIdx]
    df = df.reset_index(drop=True) # 后面需要按index索引，这边需要重置
    print(df.shape)

    #根据指定的ABBA，BABA位点数量的点计算权重
    weight = df.loc[:, ABBAname+BABAname].sum(axis=1)
    L = df.shape[0] # block数量
    num = df[ABBAname].sum(axis=1).values - df[BABAname].sum(axis=1).values
    den =  df[ABBAname].sum(axis=1).values + df[BABAname].sum(axis=1).values
    totAbba = df[ABBAname].sum(axis=1).sum()
    totBaba = df[BABAname].sum(axis=1).sum()
    totBbaa = df[BBAAname].sum(axis=1).sum()

    #block jack knife estimation
    weight = weight / np.sum(weight) #block weights
    L = len(num)
    thetaN = np.sum(num) / np.sum(den) #D statistics calculated without jackknife
    pseudo = 1 / weight
    thetaJStar = [] # D statistic in each block
    for i in range(L):
        JKD = df.loc[df.index!=i, 'Numer'].sum() / df.loc[df.index!=i, 'Denom'].sum()
        thetaJStar.append(JKD)
    meanJack = np.mean(thetaJStar)# average D statistics over the blocks
    thetaJack = L*thetaN - sum((1-weight) * thetaJStar) # jackknife D statistic
    thetaTilde = pseudo*thetaN-(pseudo-1) * thetaJStar #intermediate quantity for the variance of the jackknife D statistic
    varJack = 1/L * sum( (1/(pseudo-1)) * ((thetaTilde - thetaJack)*((thetaTilde - thetaJack))) ) #variance of the jackknife D
    Z = thetaN / np.sqrt(varJack)
    print(thetaN, thetaJack, varJack, Z, totAbba, totBaba, totBbaa, L)
    return thetaN, thetaJack, varJack, Z, totAbba, totBaba, totBbaa, L, thetaJStar




# =============================================================================
# ABBA<-c("0110","0220","0330","1001","1221","1331","2002","2112","2332","3003","3113","3223")
# BABA<-c("0101","0202","0303","1010","1212","1313","2020","2121","2323","3030","3131","3232")
# ABBAtr<-c("0110","0330","1001","1221","2112","2332","3003","3223")
# BABAtr<-c("0101","0303","1010","1212","2121","2323","3030","3232")
# BBAA<-c("0011","0022","0033","1100","1122","1133","2200","2211","2233","3300","3311","3322")
#
# =============================================================================

@click.command()
@click.option('--abbababa2', help='angsd *.abababa2')
@click.option('--reportfile', help='angsd report.Observed.txt, estAvgError.R')
@click.option('--h1')
@click.option('--h2')
@click.option('--h3')
@click.option('--h4')
@click.option('--outprefix')
def main(abbababa2, reportfile, h1, h2, h3, h4, outprefix):
    df1 = pd.read_csv(reportfile, sep='\t')
    query_index = df1.loc[(df1.H1==h1) & (df1.H2==h2) & (df1.H3==h3) & (df1.H4==h4), :].index[0]
    df = []
    with open(abbababa2) as f:
        header = f.readline().strip().split()
        for line, index in zip(f, cycle(range(df1.shape[0]))):
            if index == query_index:
                tline = [float(x) for x in line.strip().split()]
                df.append(tline)
            index += 1
    df = pd.DataFrame(df, columns=header)
    print(df.shape)
    ABBAtr = ["0110","0330","1001","1221","2112","2332","3003","3223"]
    BABAtr = ["0101","0303","1010","1212","2121","2323","3030","3232"]
    ABBA = ["0110","0220","0330","1001","1221","1331","2002","2112","2332","3003","3113","3223"]
    BABA = ["0101","0202","0303","1010","1212","1313","2020","2121","2323","3030","3131","3232"]
    BBAA = ["0011","0022","0033","1100","1122","1133","2200","2211","2233","3300","3311","3322"]
    # Observed
    print('Observed')
    thetaN, thetaJack, varJack, Z, totAbba, totBaba, totBbaa, L, thetaJStar = getJackKnife(
            df, ABBA, BABA, BBAA)
    with open(f'{outprefix}.Observed.stat', 'w') as f:
        f.write('D\tJK-D\tV(JK-D)\tZ\tnABBA\tnBABA\tnBlocks\tH1\tH2\tH3\tH4\n')
        f.write(f'{thetaN}\t{thetaJack}\t{varJack}\t{Z}\t{totAbba}\t{totBaba}\t{L}\t{h1}\t{h2}\t{h3}\t{h4}\n')
    with open(f'{outprefix}.Observed.JKDs', 'w') as f:
        for JKD in thetaJStar:
            f.write(f'{JKD}\n')

    # rmTransi
    print('rmTransi')
    thetaN, thetaJack, varJack, Z, totAbba, totBaba, totBbaa, L, thetaJStar = getJackKnife(
            df, ABBAtr, BABAtr, BBAA)
    with open(f'{outprefix}.TransRem.stat', 'w') as f:
        f.write('D\tJK-D\tV(JK-D)\tZ\tnABBA\tnBABA\tnBlocks\tH1\tH2\tH3\tH4\n')
        f.write(f'{thetaN}\t{thetaJack}\t{varJack}\t{Z}\t{totAbba}\t{totBaba}\t{L}\t{h1}\t{h2}\t{h3}\t{h4}\n')
    with open(f'{outprefix}.TransRem.JKDs', 'w') as f:
        for JKD in thetaJStar:
            f.write(f'{JKD}\n')


if __name__ == '__main__':
    main()

