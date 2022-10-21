# -*- coding: utf-8 -*-
"""
Created on Thu Aug 19 22:16:27 2021

@author: Yudongcai

@Email: yudong_cai@163.com
"""


import typer
import pysam
import subprocess
import numpy as np
import pandas as pd
from collections import defaultdict



def load_reads_info(infile, reffile=None):
    reads2AS = defaultdict(int)
    if infile.endswith('bam'):
        alnfile = pysam.AlignmentFile(infile, 'rb')
    elif infile.endswith('.cram'):
        alnfile = pysam.AlignmentFile(infile, 'rc', reference_filename=reffile)
    for seg in alnfile:
        # 不比较supplementary和secondary的情况，只去比较主要的比对结果。因为有时候把supplementary的加起来分数会高于300（两个完整的150 PE比对）
        if not (seg.is_supplementary or seg.is_secondary):
            reads2AS[seg.query_name] += seg.get_tag('AS')
    maxvalue = max(reads2AS.values())
    if maxvalue <= np.iinfo(np.uint16).max:
        reads2AS = pd.Series(reads2AS, dtype=np.uint16)
        print(f'dtype: np.uint16')
    elif maxvalue <= np.iinfo(np.uint32).max:
        reads2AS = pd.Series(reads2AS, dtype=np.uint32)
        print(f'dtype: np.uint32')
    else:
        reads2AS = pd.Series(reads2AS, dtype=np.uint64)
        print(f'dtype: np.uint64')
    return reads2AS



def selection_reads(reads2AS1, reads2AS2):
    """
    输出bam1和bam2中需要保留的reads
    """
    df = pd.concat([reads2AS1, reads2AS2], axis=1, keys=['bam1', 'bam2'], join='outer').fillna(0)

    # 保留得分大于等于在另一个BAM中的reads
    keepped_reads_bam1 = df.loc[df['bam1']>=df['bam2'], :].index
    keepped_reads_bam2 = df.loc[df['bam2']>=df['bam1'], :].index

    print(f'bam1: {reads2AS1.shape[0]:,}')
    print(f'bam1 filtered: {len(keepped_reads_bam1):,}')
    print(f'bam2: {reads2AS2.shape[0]:,}')
    print(f'bam2 filtered: {len(keepped_reads_bam2):,}')
    return keepped_reads_bam1, keepped_reads_bam2


def filter_bam(inbam, outbam, reads_set: set, reffile=None):
    outreadsIDfile = outbam + '.readsID'
    with open(outreadsIDfile, 'w') as f:
        for read in reads_set:
            f.write(f'{read}\n')
    if inbam.endswith('bam'):
        cmd_filterBAM = ['samtools', 'view',
                         '-N', outreadsIDfile,
                         '-o', outbam,
                         inbam]
    else:
        cmd_filterBAM = ['samtools', 'view',
                         '-T', reffile,
                         '-N', outreadsIDfile,
                         '-o', outbam,
                         inbam]
    runout = subprocess.run(args=cmd_filterBAM, capture_output=True)
    if runout.returncode != 0:
        print(f'filter error: {outbam}')
        print(runout.stderr)
    else:
        print(f'filtered: {outbam}')
        cmd_index = ['samtools', 'index', outbam]
        subprocess.run(args=cmd_index)
        cmd_rm = ['rm', outreadsIDfile]
        subprocess.run(args=cmd_rm)



def main(inbam1: str = typer.Argument(..., help="BAM from ref1"),
         inbam2: str = typer.Argument(..., help="BAM from ref2"),
         outbam1: str = typer.Argument(..., help="filtered BAM from ref1"),
         outbam2: str = typer.Argument(..., help="filtered BAM from ref2"),
         ref1: str = typer.Option(None, help="path of ref1, only needed for cram input"),
         ref2: str = typer.Option(None, help="path of ref2, only needed for cram input")):
    """
    bam1和bam2为同样的reads集合比对到不同的参考基因组上。
    该脚本比较每个reads在bam1和bam2中的比对分数和段数，为每一个read判断最佳的归宿，之后筛选BAM文件。
    reads在两个BAM文件中情况一样时，会在两个文件中都保留下来，也就是只对比较效果有差异的reads做分配。
    如果AS得分超过65535会有问题。三代reads的话改一下load_reads_info里面的dtype。
    """
    reads2AS1 = load_reads_info(inbam1, ref1)
    reads2AS2 = load_reads_info(inbam2, ref2)
    keepped_reads_bam1, keepped_reads_bam2 = selection_reads(reads2AS1, reads2AS2)
    del(reads2AS1)
    del(reads2AS2)
    filter_bam(inbam1, outbam1, keepped_reads_bam1, ref1)
    del(keepped_reads_bam1)
    filter_bam(inbam2, outbam2, keepped_reads_bam2, ref2)
    
    

if __name__ == '__main__':
    typer.run(main)
