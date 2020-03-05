# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 17:18:10 2019

@author: YudongCai

@Email: yudongcai216@gmail.com
"""

import click



@click.command()
@click.option('--fastq', help='fastq file')
@click.option('--outprefix', help='输出统计结果的文件前缀')
def main(fastq, outprefix):
    """
    统计之前需要先试用gff2sqliteDB.py把gff3转换为database文件
    """
    db = gffutils.FeatureDB(dbfile, keep_order=True)
    seqs = Fasta(fafile)
    gene_stat(db, seqs, outprefix)
    mRNA_stat(db, seqs, outprefix)



if __name__ == '__main__':
    main()
