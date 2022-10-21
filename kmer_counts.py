# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 15:29:54 2019

@author: YudongCai

@Email: yudongcai216@gmail.com
"""

import click
import khmer
import pandas as pd
from itertools import product


@click.command()
@click.option('--infile', help='input fasta or fastq file')
@click.option('--k', help='[k]mer', type=int)
@click.option('--outfile', help='output file')
def main(infile, k, outfile):
    kmerlist = []
    for kmer in product('ACGT', repeat=k):
        kmerlist.append(''.join(kmer))

    df = []
    seqids = []
    for seq in khmer.ReadParser(infile):
        seqids.append(seq.name)
        df.append([])
        counts = khmer.Counttable(k, 1e6, 4)
        counts.set_use_bigcount(True)
        counts.consume(seq.sequence.upper())
        for kmer in kmerlist:
            df[-1].append(counts.get(kmer))

    df = pd.DataFrame(df, columns=kmerlist, index=seqids)
    df.T.reset_index().to_csv(outfile, index=False, sep='\t')


if __name__ == '__main__':
    main()
