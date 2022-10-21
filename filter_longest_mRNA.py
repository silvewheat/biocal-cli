# -*- coding: utf-8 -*-
"""
Created on Tue Nov  6 21:20:37 2018

@author: YudongCai

@Email: yudongcai216@gmail.com
"""

from collections import defaultdict
import click


@click.command()
@click.argument('ingff')
@click.argument('outgff')
def main(ingff, outgff):
    gene2mRNA = defaultdict(list)
    mRNA2len = defaultdict(int)
    mRNA2record = defaultdict(list)
    genelist = []
    geneset = set()
    with open(ingff) as f:
        for line in f:
            tline = line.strip().split('\t')
            gtype = tline[2]
            ann = {x.split('=')[0]: x.split('=')[1] for x in tline[8].split(';')}
            if gtype == 'mRNA':
                gene_id = ann['Parent']
                if gene_id not in geneset:
                    genelist.append(gene_id)
                    geneset.add(gene_id)
                RNA_id = ann['ID']
                gene2mRNA[gene_id].append(RNA_id)
                mRNA2record[RNA_id].append(line)
            elif gtype == 'CDS':
                RNA_id = ann['Parent']
                length = int(tline[4]) - int(tline[3]) + 1
                mRNA2len[RNA_id] += length
                mRNA2record[RNA_id].append(line)
            else:
                pass
    with open(outgff, 'w') as f:
        for gene in genelist:
            mRNAs = gene2mRNA[gene]
            lengths = [mRNA2len[x] for x in mRNAs]
            maxmRNA = mRNAs[lengths.index(max(lengths))]
            f.write(''.join(mRNA2record[maxmRNA]))

if __name__ == '__main__':
    main()

