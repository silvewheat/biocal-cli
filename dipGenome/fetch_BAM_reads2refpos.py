import gzip
import typer
import pysam



def main(inbam: str = typer.Argument(..., help="BAM or cram file"),
         outfile: str = typer.Argument(..., help="output file for reads info (XX.reads2refpos.tsv.gz)"),
         reffile: str = typer.Option(None, help="path of referefile file, only needed for cram input")):
    """
    从BAM中读取reads比对到参考上的坐标位置并输出为tsv文件(gzip压缩)
    """
    with gzip.open(outfile, 'wb') as f:
        if inbam.endswith('bam'):
            alnfile = pysam.AlignmentFile(inbam, 'rb')
        elif inbam.endswith('.cram'):
            alnfile = pysam.AlignmentFile(inbam, 'rc', reference_filename=reffile)

        for seg in alnfile:
            if not (seg.is_unmapped or seg.is_duplicate or seg.is_supplementary or seg.is_qcfail or seg.is_secondary):
                # supplementary中会有hardclip导致reads的坐标发生变化
                nread = 1 if seg.is_read1 else 2
                aligned_pairs = seg.get_aligned_pairs(matches_only=True)
                querypos = ','.join([str(x[0]) for x in aligned_pairs])
                refpos = ','.join([str(x[1]) for x in aligned_pairs])
                f.write(f'{seg.query_name}\t{nread}\t{seg.reference_name}\t{querypos}\t{refpos}\n'.encode())


if __name__ == '__main__':
    typer.run(main)