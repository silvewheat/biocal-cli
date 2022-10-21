# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 11:00:15 2020

@author: YudongCai

@Email: yudongcai216@gmail.com
"""

import os
import glob
import click



@click.command()
@click.option('--suffix1')
@click.option('--suffix2')
@click.option('--path1', help='suffix1文件所在目录，默认当前目录', default='.')
@click.option('--path2', help='suffix2文件所在目录，默认当前目录', default='.')
def main(suffix1, suffix2, path1, path2):
    """
    比较当前目录下，suffix1比suffix2多出来的个体
    比如suffix1为.g.vcf.gz
    suffix2为.g.vcf.gz.tbi
    运行后可知有哪些个体的GVCF没跑完
    """
    len1 = len(suffix1) * -1
    len2 = len(suffix2) * -1
    set1 = {os.path.basename(i)[:len1] for i in glob.glob(f'{path1}/*{suffix1}')}
    set2 = {os.path.basename(i)[:len2] for i in glob.glob(f'{path2}/*{suffix2}')}
    print(f'number of suffix1: {len(set1)}')
    print(f'number of suffix2: {len(set2)}')
    for i in set1 - set2:
        print(i)


if __name__ == '__main__':
    main()
