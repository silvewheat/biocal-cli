# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 17:05:14 2019

@author: YudongCai

@Email: yudongcai216@gmail.com
"""

import click
import numpy as np
from Bio import AlignIO
from collections import Counter




@click.command()
@click.option('--maffile')
@click.option('--outfile')
def main(maffile, outfile):
    """
    输出maf比对文件
    输出每个位点的major allele类型（用参考序列（排第一个的物种）的坐标）
    MAF http://genome.ucsc.edu/FAQ/FAQformat.html#format5
    input start is zero-based
    output pos is one-based
    """
    with open(outfile, 'w') as f:
        for multiple_alignment in AlignIO.parse(maffile, "maf"):
            starts = []
            seqs = []
            for seqrec in multiple_alignment:
                starts.append((seqrec.annotations["start"]))
                seqs.append(list(str(seqrec.seq.upper())))
            seqs = np.array(seqs)
            loc = starts[0]
            for seq in seqs.T:
                if seq[0] != '-':
                    major = Counter(seq).most_common()[0][0]
                    f.write(f'{loc + 1}\t{major}\n')
                    loc += 1

if __name__ == '__main__':
    main()
