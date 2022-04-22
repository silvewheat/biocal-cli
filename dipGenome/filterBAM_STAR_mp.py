# -*- coding: utf-8 -*-
"""
Created on 2021-10-18

@author: Yudongcai

@Email: yudong_cai@163.com
"""

import os
import typer
import gzip
import time
import pysam
import subprocess
import numpy as np
import pandas as pd
from collections import defaultdict
from multiprocessing import Process, Pipe



def stats_reads_AS_nM_pipe(infile, reffile, child_conn):
    print(f'{time.asctime()} stats_reads_AS_nM: {os.path.basename(infile)}')
    reads2AS = defaultdict(int)
    reads2nM = defaultdict(int)
    if infile.endswith('bam'):
        alnfile = pysam.AlignmentFile(infile, 'rb', threads=2)
    elif infile.endswith('.cram'):
        alnfile = pysam.AlignmentFile(infile, 'rc', reference_filename=reffile, threads=2)
    for seg in alnfile:
        # markdup应该放在区分单倍体之后 #or seg.is_duplicate 
        if seg.is_unmapped  or seg.is_duplicate or seg.is_qcfail or seg.is_supplementary or seg.is_secondary:
            continue
        readname = seg.query_name
        reads2nM[readname] = seg.get_tag('nM')
        # STAR中AS得分是pair的reads加起来的得分（bwa中是分开的）
        reads2AS[readname] = seg.get_tag('AS')

    reads2AS = pd.Series(reads2AS)
    reads2nM = pd.Series(reads2nM)
    
    df = pd.concat([reads2AS, reads2nM], axis=1, keys=['AS', 'nM'], join='outer') # 不应该有缺失值
    print(df.dtypes) # dtype类型转换了说明有缺失值
    alnfile.close()
    child_conn.send(df)

def selection_reads(mdf):
    """
    输出bam1和bam2中需要保留的reads
    """
    # 根据目标位点的错配
    nreads = mdf.shape[0]
    print(f'Total reads: {nreads}')
    
    # 填充缺失值，AS缺少的话按0填充，nM缺少的话按最大值+1填充
    mdf[['AS_bam1', 'AS_bam2']].fillna(0, inplace=True, downcast='infer')
    max_nM = mdf[['nM_bam1', 'nM_bam2']].max().max()
    mdf[['nM_bam1', 'nM_bam2']].fillna(max_nM+1, inplace=True, downcast='infer')
    
    # 根据AS
    include_AS_bam1 = mdf['AS_bam1'].values > mdf['AS_bam2'].values
    include_AS_bam2 = mdf['AS_bam2'].values > mdf['AS_bam1'].values
    
    # 根据nM
    include_nM_bam1 = mdf['nM_bam1'].values < mdf['nM_bam2'].values
    include_nM_bam2 = mdf['nM_bam2'].values < mdf['nM_bam1'].values
    
    # 开始排除
    keepped_reads_bam1 = mdf.index.values[include_AS_bam1 & include_nM_bam1]
    keepped_reads_bam2 = mdf.index.values[include_AS_bam2 & include_nM_bam2]
    print(f'{len(keepped_reads_bam1)} reads keeped in BAM1.')
    print(f'{len(keepped_reads_bam2)} reads keeped in BAM2.')
    return keepped_reads_bam1, keepped_reads_bam2


def filter_bam(inbam, outbam, reads_set: set, reffile=None):
    outreadsIDfile = outbam + '.readsID'
    with open(outreadsIDfile, 'w') as f:
        for read in reads_set:
            f.write(f'{read}\n')
    # -F 1796 (UNMAP,SECONDARY,QCFAIL,DUP)
    if reffile:
        cmd_filterBAM = ['samtools', 'view',
                         '-T', reffile,
                         '-N', outreadsIDfile,
                         '-F', '1796',
                         '-o', outbam,
                         '--threads', '6',
                         '--write-index',
                         inbam]
    else:
        cmd_filterBAM = ['samtools', 'view',
                         '-N', outreadsIDfile,
                         '-F', '1796',
                         '-o', outbam,
                         '--threads', '6',
                         '--write-index',
                         inbam]
    runout = subprocess.run(args=cmd_filterBAM, capture_output=True)
    if runout.returncode != 0:
        print(f'filter error: {outbam}')
        print(runout.stderr)
    else:
        print(f'filtered: {outbam}')
        cmd_rm = ['rm', outreadsIDfile]
        subprocess.run(args=cmd_rm)


def main(inbam1: str = typer.Option(..., help="BAM from ref1"),
         inbam2: str = typer.Option(..., help="BAM from ref2"),
         outbam1: str = typer.Option(..., help="filtered BAM from ref1"),
         outbam2: str = typer.Option(..., help="filtered BAM from ref2"),
         ref1: str = typer.Option(None, help="path of ref1"),
         ref2: str = typer.Option(None, help="path of ref2"),
         saveinfo: str = typer.Option(None, help="save reads infor into this file (.parquet.gzip)")
         ):
    """
    bam1和bam2为同样的reads集合用STAR比对到不同的参考基因组上。
    该脚本根据nM和AStag去判断每个reads应该属于哪个基因组，之后筛选BAM文件。
    在两套基因组比对中没有明显差异的reads在输出结果中不保留。
    """
    parent_conn1, child_conn1 = Pipe(duplex=False)
    parent_conn2, child_conn2 = Pipe(duplex=False)
    p1 = Process(target=stats_reads_AS_nM_pipe, args=(inbam1, ref1, child_conn1))
    p2 = Process(target=stats_reads_AS_nM_pipe, args=(inbam2, ref2, child_conn2))
    p1.start()
    p2.start()
    df1 = parent_conn1.recv()
    df2 = parent_conn2.recv()
    p1.join()
    p2.join()
    mdf = pd.merge(df1, df2, left_index=True, right_index=True, how='outer', suffixes=['_bam1', '_bam2'])
    del(df1); del(df2)
    print(mdf.dtypes)
    
    keepped_reads_bam1, keepped_reads_bam2 = selection_reads(mdf)

    p1 = Process(target=filter_bam, args=(inbam1, outbam1, keepped_reads_bam1, ref1))
    p2 = Process(target=filter_bam, args=(inbam2, outbam2, keepped_reads_bam2, ref2))
    p1.start()
    p2.start()
    p1.join()
    p2.join()

    if saveinfo:
        print(f'save reads info to {saveinfo}')
        mdf.to_parquet(saveinfo, compression='gzip')  

if __name__ == '__main__':
    typer.run(main)