# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 11:05:18 2019

@author: YudongCai

@Email: yudongcai216@gmail.com
"""

import click
import gffutils
from pyfaidx import Fasta
import pandas as pd
from collections import Counter





def cal_intro(exons):
    """
    exons至少要有两个
    """
    last_end = exons[0][1]
    for exon in exons[1: ]:
        intro_len = exon[0] - last_end - 1
        last_end = exon[1]
        yield intro_len


def gene_stat(db, seqs, outprefix):
    print('start gene stats')
    types = []
    results = []
    for gene in db.features_of_type('gene'):
        chrom = gene.chrom
        start = gene.start
        end = gene.end
        GC = seqs[chrom][start-1: end].gc
        ID = gene.attributes['ID'][0]
        Name = gene.attributes['Name'][0]
        Biotype = gene.attributes['gene_biotype'][0]
        CDS = 0
        exon = 0
        mRNA = 0
        ncRNA = 0
        featuretypes = []
        exons = []
        uniq_exon = 0
        CDSs = []
        uniq_CDS = 0
        for rec in db.children(gene):
            featuretypes.append(rec.featuretype)
            types.append(rec.featuretype)
            if rec.featuretype == 'exon':
                exons.append(f'{rec.start}_{rec.end}')
            if rec.featuretype == 'CDS':
                CDSs.append(f'{rec.start}_{rec.end}')
        uniq_exon = len(set(exons))
        uniq_CDS = len(set(CDSs))
        c = Counter(featuretypes)
        CDS += c.get('CDS', 0)
        exon += c.get('exon', 0)
        mRNA += c.get('mRNA', 0)
        ncRNA += c.get('ncRNA', 0)
        results.append([ID, Name, Biotype, GC, mRNA, ncRNA, exon, uniq_exon, CDS, uniq_CDS])
    types = Counter(types)
    print(types)
    df = pd.DataFrame(results, columns=['ID', 'Name', 'Biotype', 'GC', 'mRNA', 'ncRNA', 'exon', 'uniq_exon', 'CDS', 'uniq_CDS'])
    df.to_csv(f'{outprefix}.gene.tsv', sep='\t', index=False, float_format='%.4f')

def mRNA_stat(db, seqs, outprefix):
    print('start mRNA stats')
    types = []
    results = []
    exon_lens = []
    intro_lens = []
    for mRNA in db.features_of_type('mRNA'):
        chrom = mRNA.chrom
        start = mRNA.start
        end = mRNA.end
        GC = seqs[chrom][start-1: end].gc
        ID = mRNA.attributes['ID'][0]
        geneName = mRNA.attributes['gene'][0]
        transcript_id = mRNA.attributes['transcript_id'][0]
        CDS = 0
        exon = 0
        featuretypes = []
        exons = []
        for rec in db.children(mRNA, order_by='start'):
            featuretypes.append(rec.featuretype)
            types.append(rec.featuretype)
            if rec.featuretype == 'exon':
                start, end = rec.start, rec.end
                exons.append([start, end])
                exon_lens.append(end - start + 1)
            if len(exons) > 1:
                for intro_len in cal_intro(exons):
                    intro_lens.append(intro_len)
        c = Counter(featuretypes)
        CDS += c.get('CDS', 0)
        exon += c.get('exon', 0)
        results.append([ID, geneName, transcript_id, GC, exon, CDS])
    types = Counter(types)
    print(types)
    df = pd.DataFrame(results, columns=['ID', 'geneName', 'transcript_id', 'GC', 'exon', 'CDS'])
    df.to_csv(f'{outprefix}.mRNA.tsv', sep='\t', index=False, float_format='%.4f')
    with open(f'{outprefix}.exon_len.tsv', 'w') as f:
        for exon_len in exon_lens:
            f.write(f'{exon_len}\n')
    with open(f'{outprefix}.intro_len.tsv', 'w') as f:
        for intro_len in intro_lens:
            f.write(f'{intro_len}\n')



@click.command()
@click.option('--dbfile', help='用gffutils生成的sqlite database')
@click.option('--fafile', help='gff对应的fasta文件')
@click.option('--outprefix', help='输出统计结果的文件前缀')
def main(dbfile, fafile, outprefix):
    """
    统计之前需要先试用gff2sqliteDB.py把gff3转换为database文件
    """
    db = gffutils.FeatureDB(dbfile, keep_order=True)
    seqs = Fasta(fafile)
    gene_stat(db, seqs, outprefix)
    mRNA_stat(db, seqs, outprefix)



if __name__ == '__main__':
    main()




