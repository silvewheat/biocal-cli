# -*- coding: utf-8 -*-
"""
Created on 2022 10-14

@author: Yudongcai

@Email: yudong_cai@163.com
"""
import gzip
import typer
import numpy as np
import pandas as pd
from collections import defaultdict

def load_fa(fafile):
    seqdict = {}
    with open(fafile) as f:
        seq_id = f.readline().strip().split()[0][1:]
        tmp_seq = []
        for line in f:
            if line[0] != '>':
                tmp_seq.append(line.strip())
            else:
                seqdict[seq_id] = ''.join(tmp_seq)
                seq_id = line.strip().split()[0][1:]
                tmp_seq = []
        else:
            seqdict[seq_id] = ''.join(tmp_seq)
    return seqdict

def load_anc(ancfile):
    aa = defaultdict(dict)
    assert ancfile[-3:] == '.gz'
    with gzip.open(ancfile, 'rb') as f:
        for line in f:
            chrom, pos, anc_allel = line.decode().strip().split()
            aa[chrom][int(pos)-1] = anc_allel # convert to 0-based
    return aa

def replace_fa(seq, aa, outfasta):
    chroms = [x for x in seq.keys() if x in aa.keys()]
    fa_length = 0
    num_anc_info = 0
    num_replaced = 0
    with open(outfasta, 'w') as f:
        for ctg in chroms:
            print(ctg)
            site2allele = aa[ctg]
            content = list(seq[ctg])
            fa_length += len(content)
            num_anc_info += len(site2allele)
            for pos, anc_allel in site2allele.items():
                if content[pos].upper() == anc_allel.upper():
                    pass
                else:
                    if content[pos] != 'N':
                        content[pos] = anc_allel
                        num_replaced += 1
            f.write(f'>{ctg}\n')
            f.write(''.join(content) + '\n')
    print(f'fa_length: {fa_length}\nnum_anc_info: {num_anc_info}\nnum_replaced: {num_replaced}')

def main(infasta,
         ancfile,
         outfasta):
    seq = load_fa(infasta)
    aa = load_anc(ancfile)
    replace_fa(seq, aa, outfasta)

if __name__ == '__main__':
    typer.run(main)