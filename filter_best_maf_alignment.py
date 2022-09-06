# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 15:21:02 2021

@author: Yudongcai

@Email: yudong_cai@163.com
"""

import re
import typer
import numpy as np
import pandas as pd


def load_maf(infile, max_mismap):
    re_score = re.compile(r'score=(\S+)\b')
    re_mismap = re.compile(r'mismap=(\S+)\b')
    df = []
    n_seg = 0
    with open(infile) as f:
        for line in f:
            if line[0] == 'a':
                n_seg += 1
                score = int(re_score.search(line).group(1))
                mismap = float(re_mismap.search(line).group(1))
                if mismap > max_mismap:
                    continue
                ref_line = f.readline()
                line_tag, chrom, start, aln_length, strand = ref_line.strip().split()[:5]
                if strand == '+':
                    df.append([n_seg, chrom, int(start), int(start)+int(aln_length), score, mismap])
                elif strand == '-':
                    df.append([n_seg, chrom, int(start)-int(aln_length), int(start), score, mismap])
                else:
                    raise Exception
    df = pd.DataFrame(df, columns=['n_seg', 'chrom', 'start', 'end', 'score', 'mismap'])
    df['mismap_neg_log10'] = np.log10(df['mismap'].values)*-1
    df = df.sort_values(['mismap_neg_log10', 'score'], ascending=False)
    return df


def filter_aln_seg(df, min_overhang_rate):
    n_seg_keep = []
    chrom_len = df['end'].max()
    seq = np.zeros(chrom_len, dtype=np.float16)
    for seg in df[['start', 'end', 'mismap_neg_log10', 'n_seg']].itertuples():
        # 判断存在没有比对覆盖的区域
        if seq[seg.start: seg.end].min() == 0:
            seg_length = seg.end - seg.start
            values1 = seq[seg.start: seg.end]
            # 判断与以有比对区域重叠部分的mismap，不能比已有的mismap值低
            # 判断超出已有比对部分的区域，不能低于当前比对片段的一定比例
            overlaped_length = np.sum(values1 > 0)
            overhang_rate = (seg_length - overlaped_length) / seg_length
            if not overhang_rate >= min_overhang_rate:
                continue
#            overlaped_mean_values = 0 if overlaped_length == 0 else np.mean(values1[values1 > 0])
            overlaped_max_values = values1.max()
#            if not seg.mismap_neg_log10 >= overlaped_mean_values:
            if not seg.mismap_neg_log10 >= overlaped_max_values:
                continue
            values2 = np.full(values1.shape, seg.mismap_neg_log10)
            seq[seg.start: seg.end] = np.max([values1, values2], axis=0) # 只填充新的部分
            n_seg_keep.append(seg.n_seg)
    return set(n_seg_keep)


def filter_maf_file(infile, outfile, n_seg_keep):
    n_seg = 0
    with open(infile) as f1, open(outfile, 'w') as f2:
        for line in f1:
            if line[0] == 'a':
                n_seg += 1
            if n_seg in n_seg_keep:
                f2.write(line)


def main(infile: str = typer.Argument(..., help="input MAF file, reference can only contain one chromsome"),
         outfile: str = typer.Argument(..., help="output filtered MAF file"),
         max_mismap: float = typer.Option(1e-5, help="mismap cutoff, only keep segments lower than this value"),
         min_overhang_rate: float = typer.Option(0.25, help="超出已有比对片段覆盖的区域至少一定比例")):
    """
    从LAST比对结果中提取最佳比对
    贪婪算法的思想，优先选取高质量比对铺满整个target/reference序列
    """
    df = load_maf(infile, max_mismap)
    assert df['chrom'].unique().shape[0] == 1 # maf文件只能包含一条reference的contig
    n_seg_keep = filter_aln_seg(df, min_overhang_rate)
    filter_maf_file(infile, outfile, n_seg_keep)


if __name__ == '__main__':
    typer.run(main)