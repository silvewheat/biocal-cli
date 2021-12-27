# -*- coding: utf-8 -*-
"""
Created on 2021-09-08

@author: Yudongcai

@Email: yudong_cai@163.com
"""

import typer
import gzip
import pyfastx
import mappy as mp
import numpy as np



def main(insites: str = typer.Argument(..., help="带转换的来自queryfa的位点，gzip压缩，两列chrom\tpos，可以是vcf.gz"),
         queryfa: str = typer.Argument(..., help="query fasta file, insites对应的fasta"),
         targetfa: str = typer.Argument(..., help="target fasta file, 目标fasta"),
         outfile: str = typer.Argument(..., help='output file (.tsv.gz)')):
    """把insites里的位点从queryfa对应到targetfa，基于minimap2比对"""
    a = mp.Aligner(targetfa, preset='sr')
    if not a: raise Exception("ERROR: failed to load/build index")
    
    fa = pyfastx.Fasta(queryfa)
    
    with gzip.open(insites, 'rb') as f1, gzip.open(outfile, 'wb') as f2:
        f2.write('qchrom\tqpos\ttchrom\ttpos\n'.encode())
        for line in f1:
            tline = line.decode().split()
            if tline[0][0] != '#':
                chrom, pos = tline[:2]
                pos = int(pos)
                seq1 = fa.fetch(chrom, (max(0, pos-150), pos-1)) # 开头几个pos可能小于150,就用0替代
                seq2 = fa.fetch(chrom, (pos+1, pos+150))
                
                locs = []
                ctgs = []
                for hit in a.map(seq1, seq2):
                    if hit.is_primary:
                        locs.append(hit.r_st)
                        locs.append(hit.r_en)
                        ctgs.append(hit.ctg)
                if len(ctgs)==2 and (len(set(ctgs))==1):
                    locs.sort()
                    if locs[2] - locs[1] == 1:
                        f2.write(f'{chrom}\t{pos}\t{ctgs[0]}\t{locs[1] + 1}\n'.encode())

if __name__ == '__main__':
    typer.run(main)
    