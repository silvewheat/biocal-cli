# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 09:17:27 2019

@author: YudongCai

@Email: yudongcai216@gmail.com
"""

import click
import numpy as np
from Bio import AlignIO
from collections import Counter




@click.command()
@click.option('--maffile')
@click.option('--query', help='query基因组的名字')
@click.option('--outfile')
def main(maffile, outfile):
    """
    输出maf比对文件
    输出每个位点的major allele类型（只包含在ref和query基因组中都存在的位点）
    MAF http://genome.ucsc.edu/FAQ/FAQformat.html#format5
    input start is zero-based
    output pos is one-based
    maf中染色体命名格式为 基因组名.染色体号
    """
query_genome = 'panTro2'
with open(outfile, 'w') as f:
    f.write('refGenome\trefChrom\trefLoc\trefStrand\trefNuc\tqueryGenome\tqueryChom\tqueryLoc\tqueryStrand\tqueryNuc\tMajorAllele\n')
    for multiple_alignment in AlignIO.parse(maffile, "maf"):
        starts = []
        seqs = []
        genomes = []
        strands = []
        chroms = []
        for seqrec in multiple_alignment:
            starts.append((seqrec.annotations["start"]))
            seqs.append(list(str(seqrec.seq.upper())))
            genomes.append(seqrec.id.split('.')[0]) # 物种.染色体
            chroms.append(seqrec.id.split('.')[-1])
            strands.append(seqrec.annotations['strand'])
        if query_genome in genomes:
            seqs = np.array(seqs)
            ref_loc = starts[0]
            print(genomes)
            query_index = genomes.index(query_genome)
            query_loc = starts[query_index]
            for seq in seqs.T:
                if seq[0] != '-':
                    major = Counter(seq).most_common()[0][0]
                    ref_loc += strands[0] # 1正链，-1负链
                    if seq[query_index] != '-':
                        query_loc += strands[query_index]
                        f.write(f'{genomes[0]}\t{chroms[0]}\t{ref_loc + 1}\t{strands[0]}\t{seq[0]}\t{genomes[query_index]}\t{chroms[query_index]}\t{query_loc+1}\t{strands[query_index]}\t{seq[query_index]}\t{major}\n')
if __name__ == '__main__':
    main()
