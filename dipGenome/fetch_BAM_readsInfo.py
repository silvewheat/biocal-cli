# -*- coding: utf-8 -*-
"""
Created on 2021-09-07

@author: Yudongcai

@Email: yudong_cai@163.com
"""

import gzip
import typer
import pysam



def main(inbam: str = typer.Argument(..., help="BAM or cram file"),
         outfile: str = typer.Argument(..., help="output file for reads info (XX.readsInfo.tsv.gz)"),
         reffile: str = typer.Option(None, help="path of referefile file, only needed for cram input")):
    """
    读取reads比对信息并输出为tsv文件(gzip压缩)
    """
    with gzip.open(outfile, 'wb') as f:
#        f.write('header'.encode())
        if inbam.endswith('bam'):
            alnfile = pysam.AlignmentFile(inbam, 'rb')
        elif inbam.endswith('.cram'):
            alnfile = pysam.AlignmentFile(inbam, 'rc', reference_filename=reffile)
        for seg in alnfile:
            readname = seg.query_name
            is_sup = int(seg.is_supplementary)
            is_sec = int(seg.is_secondary)
            is_umap = int(seg.is_unmapped)
            if seg.is_read1:
                read12 = 'R1'
            elif seg.is_read2:
                read12 = 'R2'
            else:
                read12 = ''
            
            if seg.is_supplementary:
                supOrSec = 'Sup'
            elif seg.is_secondary:
                supOrSec = 'Sec'
            else:
                supOrSec = ''

            NM = seg.get_tag('NM')
            AS = seg.get_tag('AS')
            f.write(f'{readname}\t{read12}\t{supOrSec}\t{NM}\t{AS}\n'.encode())

if __name__ == '__main__':
    typer.run(main)