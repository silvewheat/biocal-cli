# -*- coding: utf-8 -*-
"""
Created on 2021-09-14

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




def timer_func(func):
    # This function shows the execution time of
    # the function object passed
    def wrap_func(*args, **kwargs):
        t1 = time.time()
        result = func(*args, **kwargs)
        t2 = time.time()
        runtime = t2 - t1
        h = runtime // 3600
        m = (runtime - h*3600) // 60
        s = runtime - h*3600 - m*60
        print(f'Function {func.__name__!r} executed in {h:.0f}h {m:.0f}m {s:.0f}s ({runtime}s)')
        return result
    return wrap_func

def select_min_dtype_uint(value: int):
    """选择最小的uint类型，无负数"""
    assert value >= 0
    if value <= np.iinfo(np.uint8).max:
        return np.uint8
    elif value <= np.iinfo(np.uint16).max:
        return np.uint16
    elif value <= np.iinfo(np.uint32).max:
        return np.uint32
    elif value <= np.iinfo(np.uint64).max:
        return np.uint64
    else:
        return int

@timer_func
def stat_query_mismatch_pipe(alnfile: str, reffile: str, sitesfile: str, child_conn):
    """
    统计指定位点的错配情况
    alnfile为bam/cram文件
    reffile为对应的参考基因组文件
    sitesfile为待检测的位点，三列，gzip或bgzip压缩，第一列为染色体，第二列为坐标位置（1-based），第三列为对应染色体上的碱基
    """
    print(f'{time.asctime()} stat_query_mismatch_pipe: {os.path.basename(alnfile)}, {os.path.basename(sitesfile)}')

    # 读取位点信息
    chrom2sites = defaultdict(list)
    chrom2sites2refbase = defaultdict(dict)
    with gzip.open(sitesfile, 'rb') as f:
        for line in f:
            tline = line.decode().strip().split()
            if tline[0][0] != '#':
                chrom = tline[0]
                pos = int(tline[1])-1 # 后续pysam输入的是0-base索引，在这儿提前转换了
                refbase = tline[2]
                chrom2sites[chrom].append(pos)
                chrom2sites2refbase[chrom][pos] = refbase
    
    # 打开bam文件
    if alnfile.endswith('bam'):
        alnfile = pysam.AlignmentFile(alnfile, 'rb', threads=2)
    elif alnfile.endswith('.cram'):
        alnfile = pysam.AlignmentFile(alnfile, 'rc', reference_filename=reffile, threads=2)
    
    # 开始遍历位点进行判断
    reads2QMis = defaultdict(int) # 每对reads在目标位点(query sites)上与参考基因组(reffile)之间的错配数量。统计错配的话，后续缺失值用0填充就比较合理。
    reads2nQuery = defaultdict(int) # 记录每个reads cover到了几个目标位点，用来区分没有错配还是没有cover到query位点
    
    for chrom, sites in chrom2sites.items():
        sites_set = set(sites)
        print(f'chrom: {chrom}, nsites: {len(sites)}')
        sites2refbase = chrom2sites2refbase[chrom]
        # stepper='all': skip reads in which any of the following flags are set: BAM_FUNMAP, BAM_FSECONDARY, BAM_FQCFAIL, BAM_FDUP
        for pileupcolumn in alnfile.pileup(contig=chrom, truncate=True, stepper='all', ignore_orphans=False,
                                           min_base_quality=0, min_mapping_quality=0):
            if pileupcolumn.reference_pos in sites_set:
                refbase = sites2refbase[pileupcolumn.reference_pos].upper()
                for pileupread in pileupcolumn.pileups:
                    readname = pileupread.alignment.query_name
                    reads2nQuery[readname] += 1
                    if not pileupread.is_del:
                        # query position is None if is_del or is_refskip is set. 统一大小写进行比较，fasta前面已经转过了
                        if pileupread.alignment.query_sequence[pileupread.query_position].upper() != refbase:
                            reads2QMis[readname] += 1
                    else:
                        reads2QMis[readname] += 1
    alnfile.close()
    
    maxnQuery = max(reads2nQuery.values())
    nQuery_dtype = select_min_dtype_uint(maxnQuery)
    reads2nQuery = pd.Series(reads2nQuery, dtype=nQuery_dtype)
    print(f'{time.asctime()} max(nQuery): {maxnQuery}, select dtype: {nQuery_dtype}')
    
    maxQMis = max(reads2QMis.values())
    QMis_dtype = select_min_dtype_uint(maxQMis)
    reads2QMis = pd.Series(reads2QMis, dtype=QMis_dtype)
    print(f'{time.asctime()} max(QMis): {maxQMis}, select dtype: {QMis_dtype}')
    
    child_conn.send((reads2nQuery, reads2QMis))


@timer_func
def stats_reads_AS_NM_pipe(infile, reffile, child_conn):
    print(f'{time.asctime()} stats_reads_AS_NM: {os.path.basename(infile)}')
    reads2AS = defaultdict(int)
    reads2NM = defaultdict(int)
    if infile.endswith('bam'):
        alnfile = pysam.AlignmentFile(infile, 'rb', threads=2)
    elif infile.endswith('.cram'):
        alnfile = pysam.AlignmentFile(infile, 'rc', reference_filename=reffile, threads=2)
    for seg in alnfile:
        if seg.is_unmapped or seg.is_duplicate or seg.is_qcfail or seg.is_secondary:
            continue
        readname = seg.query_name
        reads2NM[readname] += seg.get_tag('NM')
        # 对于AS不比较supplementary和secondary的情况，只去比较主要的比对结果。因为有时候把supplementary的加起来分数会高于300（两个完整的150 PE比对）
        if not seg.is_supplementary:
            reads2AS[readname] += seg.get_tag('AS')
        else:
            reads2AS[readname] += 0 # 保证reads2NM和reads2AS的key一致
    maxAS = max(reads2AS.values())
    AS_dtype = select_min_dtype_uint(maxAS)
    reads2AS = pd.Series(reads2AS, dtype=AS_dtype)
    print(f'max(AS): {maxAS}, select dtype: {AS_dtype}')
    
    maxNM = max(reads2NM.values())
    NM_dtype = select_min_dtype_uint(maxNM)
    reads2NM = pd.Series(reads2NM, dtype=NM_dtype)
    print(f'max(NM): {maxNM}, select dtype: {NM_dtype}')
    
    df = pd.concat([reads2AS, reads2NM], axis=1, keys=['AS', 'NM'], join='outer') # 不应该有缺失值
    print(df.dtypes) # dtype类型转换了说明有缺失值
    alnfile.close()
    child_conn.send(df)

@timer_func
def selection_reads(mdf):
    """
    输出bam1和bam2中需要保留的reads
    """
    # 根据目标位点的错配
    nreads = mdf.shape[0]
    print(f'Total reads: {nreads}')
    # 因为在hap1和hap2的基因组中找到的query位点数量可能不同，做标准化也不合适，还是得分母相等再比较分子
    selection_both_cover = (mdf['nQuery_bam1'].values > 0) & (mdf['nQuery_bam2'].values > 0)
    selection_same_cover = (mdf['nQuery_bam1'].values == mdf['nQuery_bam2'].values) & selection_both_cover
    exclude_queryMisMatch_bam1 = (mdf['QMis_bam1'].values > mdf['QMis_bam2'].values) & selection_same_cover
    exclude_queryMisMatch_bam2 = (mdf['QMis_bam1'].values < mdf['QMis_bam2'].values) & selection_same_cover
#    exclude_queryMisMatch_bam1 = mdf['QMis_bam1'].values[selection_same_cover] > mdf['QMis_bam2'].values[selection_same_cover]
#    exclude_queryMisMatch_bam2 = mdf['QMis_bam1'].values[selection_same_cover] < mdf['QMis_bam2'].values[selection_same_cover]
    
#    ay_qmisRate_bam1 = np.full(nreads, np.nan)
#    selection_nQuery = mdf['nQuery_bam1'].values > 0
#    ay_qmisRate_bam1[selection_nQuery] = mdf['QMis_bam1'].values[selection_nQuery] / mdf['nQuery_bam1'].values[selection_nQuery]
#    
#    ay_qmisRate_bam2 = np.full(nreads, np.nan)
#    selection_nQuery = mdf['nQuery_bam2'].values > 0
#    ay_qmisRate_bam2[selection_nQuery] = mdf['QMis_bam2'].values[selection_nQuery] / mdf['nQuery_bam2'].values[selection_nQuery]
#    
#    exclude_queryMisMatch_bam1 = ay_qmisRate_bam1 > ay_qmisRate_bam2
#    exclude_queryMisMatch_bam2 = ay_qmisRate_bam1 < ay_qmisRate_bam2
    
    # 根据AS
    exclude_AS_bam1 = mdf['AS_bam1'].values < mdf['AS_bam2'].values
    exclude_AS_bam2 = mdf['AS_bam1'].values > mdf['AS_bam2'].values
    
    # 根据NM
    exclude_NM_bam1 = mdf['NM_bam1'].values > mdf['NM_bam2'].values
    exclude_NM_bam2 = mdf['NM_bam1'].values < mdf['NM_bam2'].values
    
    # 开始排除
    keepped_reads_bam1 = set(mdf.index)
    keepped_reads_bam1 -= set(mdf[exclude_queryMisMatch_bam1].index)
    print(f'BAM1 after exclude_queryMisMatch: {len(keepped_reads_bam1)}')
    keepped_reads_bam1 -= set(mdf[exclude_AS_bam1].index)
    print(f'BAM1 after exclude_AS: {len(keepped_reads_bam1)}')
    keepped_reads_bam1 -= set(mdf[exclude_NM_bam1].index)
    print(f'BAM1 after exclude_NM: {len(keepped_reads_bam1)}')
    
    keepped_reads_bam2 = set(mdf.index)
    keepped_reads_bam2 -= set(mdf[exclude_queryMisMatch_bam2].index)
    print(f'BAM2 after exclude_queryMisMatch: {len(keepped_reads_bam2)}')
    keepped_reads_bam2 -= set(mdf[exclude_AS_bam2].index)
    print(f'BAM2 after exclude_AS: {len(keepped_reads_bam2)}')
    keepped_reads_bam2 -= set(mdf[exclude_NM_bam2].index)
    print(f'BAM2 after exclude_NM: {len(keepped_reads_bam2)}')
    return keepped_reads_bam1, keepped_reads_bam2

@timer_func
def filter_bam(inbam, outbam, reads_set: set, reffile=None):
    outreadsIDfile = outbam + '.readsID'
    with open(outreadsIDfile, 'w') as f:
        for read in reads_set:
            f.write(f'{read}\n')
    # -F 1796 (UNMAP,SECONDARY,QCFAIL,DUP)
    cmd_filterBAM = ['samtools', 'view',
                     '-T', reffile,
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
         ref1: str = typer.Option(..., help="path of ref1"),
         ref2: str = typer.Option(..., help="path of ref2"),
         sitesfile1: str = typer.Option(..., help="het sites file for ref1. chrom\tpos, gzipped, can be vcf.gz"),
         sitesfile2: str = typer.Option(..., help="het sites file for ref2. chrom\tpos, gzipped, can be vcf.gz"),
         saveinfo: str = typer.Option(None, help="save reads infor into this file (.tsv.gz)")
         ):
    """
    bam1和bam2为同样的reads集合比对到不同的参考基因组上。
    该脚本根据已知的两个单倍体基因组之间的差异位点(sitesfile1 和 sitesfile2)去判断每个reads应该属于哪个基因组，之后筛选BAM文件。
    reads在两个BAM文件中情况一样时，会在两个文件中都保留下来，也就是只对比较效果有差异的reads做分配。
    只适用于二代reads，bwa mapping。
    """
    parent_conn1, child_conn1 = Pipe(duplex=False)
    parent_conn2, child_conn2 = Pipe(duplex=False)
    p1 = Process(target=stats_reads_AS_NM_pipe, args=(inbam1, ref1, child_conn1))
    p2 = Process(target=stats_reads_AS_NM_pipe, args=(inbam2, ref2, child_conn2))
    p1.start()
    p2.start()
    df1 = parent_conn1.recv()
    df2 = parent_conn2.recv()
    p1.join()
    p2.join()
    mdf = pd.merge(df1, df2, left_index=True, right_index=True, how='outer', suffixes=['_bam1', '_bam2'])
    del(df1); del(df2)
    print(mdf.dtypes)

    p1 = Process(target=stat_query_mismatch_pipe, args=(inbam1, ref1, sitesfile1, child_conn1))
    p2 = Process(target=stat_query_mismatch_pipe, args=(inbam2, ref2, sitesfile2, child_conn2))
    p1.start()
    p2.start()
    reads2nQuery1, reads2QMis1 = parent_conn1.recv()
    reads2nQuery2, reads2QMis2 = parent_conn2.recv()
    p1.join()
    p2.join()
    
    mdf['QMis_bam1'] = mdf.index.map(lambda x: reads2QMis1.get(x, 0)); del(reads2QMis1)
    mdf['QMis_bam2'] = mdf.index.map(lambda x: reads2QMis2.get(x, 0)); del(reads2QMis2)
    mdf['nQuery_bam1'] = mdf.index.map(lambda x: reads2nQuery1.get(x, 0)); del(reads2nQuery1)
    mdf['nQuery_bam2'] = mdf.index.map(lambda x: reads2nQuery2.get(x, 0)); del(reads2nQuery2)
    
    keepped_reads_bam1, keepped_reads_bam2 = selection_reads(mdf)

    p1 = Process(target=filter_bam, args=(inbam1, outbam1, keepped_reads_bam1, ref1))
    p2 = Process(target=filter_bam, args=(inbam2, outbam2, keepped_reads_bam2, ref2))
    p1.start()
    p2.start()
    p1.join()
    p2.join()

    if saveinfo:
        print(f'save reads info to {saveinfo}')
        mdf.to_csv(saveinfo, sep='\t')

if __name__ == '__main__':
    typer.run(main)