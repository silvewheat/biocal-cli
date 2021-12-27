import os
import time
import gzip
import typer
import pysam
import numpy as np
import pandas as pd
from scipy.sparse import lil_matrix, save_npz


def timer_func(func):
    # This function shows the execution time of
    # the function object passed
    def wrap_func(*args, **kwargs):
        t1 = time.time()
        print(f'{time.asctime()}: {func.__name__!r} start.')
        result = func(*args, **kwargs)
        t2 = time.time()
        runtime = t2 - t1
        h = runtime // 3600
        m = (runtime - h*3600) // 60
        s = runtime - h*3600 - m*60
        print(f'{time.asctime()}: {func.__name__!r} finished.')
        print(f'Function {func.__name__!r} executed in {h:.0f}h {m:.0f}m {s:.0f}s ({runtime}s)')
        return result
    return wrap_func

              
@timer_func
def load_refinfo(alnfile):
    with pysam.AlignmentFile(alnfile, 'rb') as alnfile:
        offsets = []
        rlens = np.array(alnfile.lengths)
        rnames = np.array(alnfile.references)
        chrom2length = {}
        chrom2offset = {}
        offset = 0
        for chrom, length in zip(rnames, rlens):
            chrom2length[chrom] = length
            chrom2offset[chrom] = offset
            offsets.append(offset)
            offset += length
    offsets = np.array(offsets)
    return chrom2offset, chrom2length, rnames, offsets


def find_read(name_indexed, readname, nread):
    """
    readname为read名字
    nread为1或2，表示read1还是read2
    """
    try:
        for seg in name_indexed.find(readname):
            if not (seg.is_unmapped or seg.is_duplicate or seg.is_supplementary or seg.is_qcfail or seg.is_secondary):
                nread_current = 1 if seg.is_read1 else 2
                if nread_current == nread:
                    return seg
        else:
            return None
    except KeyError:
        return None

@timer_func
def load_pos_to_matrix(alnfile_hap1, name_indexed_hap2, ctg_hap1, chrom2length_hap1, chrom2length_hap2, chrom2offset_hap2):
    """
    根据reads的比对关系，把
    """
    # matrix是n×m的lil稀疏矩阵, n为当前query染色体(ctg_hap1)长度, m为target全基因组长度
    n = chrom2length_hap1[ctg_hap1]
    m = sum(chrom2length_hap2.values())
    print(ctg_hap1, n, m)
    matrix_R2T = lil_matrix((n, m), dtype=np.uint16)
    count_reads = 0
    with pysam.AlignmentFile(alnfile_hap1, 'rb') as alnfile_hap1:
        for seg_hap1 in alnfile_hap1.fetch(contig=ctg_hap1):
            count_reads += 1
            if count_reads % 1000000 == 0:
                print(f'{count_reads//1000000:.0f} M reads processed.')
            if not (seg_hap1.is_unmapped or seg_hap1.is_duplicate or seg_hap1.is_supplementary or seg_hap1.is_qcfail or seg_hap1.is_secondary):
                # mapping quality 低于30的不要
                if seg_hap1.mapping_quality >= 30:
                    nread = 1 if seg_hap1.is_read1 else 2
                    seg_hap2 = find_read(name_indexed_hap2, seg_hap1.query_name, nread)
                    # mapping quality 低于30的不要
                    if seg_hap2 and seg_hap2.mapping_quality >= 30:
                        offset_hap2 = chrom2offset_hap2[seg_hap2.reference_name]
                        pairs = pd.concat([pd.Series(dict(seg_hap1.get_aligned_pairs(matches_only=True))), pd.Series(dict(seg_hap2.get_aligned_pairs(matches_only=True)))], axis=1).dropna().astype(int).values
                        for x, y in pairs:
                            matrix_R2T[x, offset_hap2+y] += 1
    return matrix_R2T.tocsr()

@timer_func
def matrix_to_pos(ctg_hap1, matrix, outfile, rnames_hap2, offsets_hap2, chrom2offset_hap2):
    """
    matrix: scipy.sparse.csr.csr_matrix
    """
    with gzip.open(outfile, 'ab') as f:
        pos_hap1 = -1
        for max_count, idx_hap2, sum_count in zip(matrix.max(axis=1).toarray(), matrix.argmax(axis=1), matrix.sum(axis=1)):
            pos_hap1 += 1
            max_count = max_count[0]
            idx_hap2 = idx_hap2[0, 0]
            sum_count = sum_count[0, 0]
            if sum_count > 0:
                perct = max_count / sum_count
                ctg_hap2 = rnames_hap2[(idx_hap2 - offsets_hap2) >= 0][-1]
                offset_hap2 = chrom2offset_hap2[ctg_hap2]
                pos_hap2 = idx_hap2 - offset_hap2
                
                f.write(f'{ctg_hap1}\t{pos_hap1+1}\t{ctg_hap2}\t{pos_hap2+1}\t{perct:.3f}\t{sum_count}\n'.encode())

@timer_func
def write_header(outfile):
    with gzip.open(outfile, 'wb') as f:
        f.write(f'chrom_query\tpos_query\tchrom_target\tpos_target\tperct_support\tcount_sum\n'.encode())



def main(inbam1: str = typer.Option(..., help="BAM from ref1"),
         inbam2: str = typer.Option(..., help="BAM from ref2"),
         ctg_hap1: str = typer.Option(None, help="link this contig in ref1"),
         ctgs_hap1: str = typer.Option(None, help="link these contigs in ref1 (list file)"),
         outdir: str = typer.Option(..., help="output link file from ref1 to ref2"),
         savematrix: str = typer.Option(None, help="save link matrix to this file")):
    """根据reads比对情况，输出ref1到ref2的对应关系。
    只支持BAM文件。
    """
    t1 = time.time()
    print('load reference infor')
    chrom2offset_hap1, chrom2length_hap1, rnames_hap1, offsets_hap1 = load_refinfo(inbam1)
    chrom2offset_hap2, chrom2length_hap2, rnames_hap2, offsets_hap2 = load_refinfo(inbam2)
    
    print('build reads index')
    alnfile_hap2 = pysam.AlignmentFile(inbam2, 'rb')
    name_indexed_hap2 = pysam.IndexedReads(alnfile_hap2)
    name_indexed_hap2.build()
    
    print('start link:')
    if ctgs_hap1:
        ctgs_hap1 = [x.strip() for x in open(ctgs_hap1)]
    else:
        ctgs_hap1 = [ctg_hap1, ]
    for ctg_hap1 in ctgs_hap1:
        print(f'ctg: {ctg_hap1}')
        matrix_R2T = load_pos_to_matrix(inbam1, name_indexed_hap2, ctg_hap1, chrom2length_hap1, chrom2length_hap2, chrom2offset_hap2)
        if savematrix:
            save_npz(os.path.join(savematrix, f'{ctg_hap1}.npz'), matrix_R2T)
#        matrix_to_pos(ctg_hap1, matrix_R2T, os.path.join(outdir, f'{ctg_hap1}.linkpos.tsv.gz'), rnames_hap2, offsets_hap2, chrom2offset_hap2)


    t2 = time.time()
    runtime = t2 - t1
    h = runtime // 3600
    m = (runtime - h*3600) // 60
    s = runtime - h*3600 - m*60
    print(f'Finished! Executed in {h:.0f}h {m:.0f}m {s:.0f}s ({runtime}s)')


if __name__ == '__main__':
    typer.run(main)