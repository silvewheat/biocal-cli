# -*- coding: utf-8 -*-
"""
Created on Sun Jan 27 14:42:14 2019

@author: YudongCai

@Email: yudongcai216@gmail.com
"""

import sys
import click
import pysam
import numpy as np
import pandas as pd
from collections import Counter



def load_alnfile(alnfile, reffile):
    if alnfile.split('.')[-1] == 'cram':
        if not reffile:
            sys.exit("cram file need refrence fasta in --reffile")
        aln = pysam.AlignmentFile(alnfile, reference_filename=reffile)
    elif alnfile.split('.')[-1] == 'bam':
        aln = pysam.AlignmentFile(alnfile, 'rb')
    elif alnfile.split('.')[-1] == 'sam':
        aln = pysam.AlignmentFile(alnfile, 'r')
    else:
        sys.exit("file suffix")
    return aln



def load_sites(sitesfile):
    """
    四列，没header
    contig pos ref alt
    """
    df = pd.read_csv(sitesfile, sep='\s+', header=None, names=['contig', 'pos', 'ref', 'alt'])
    return df



@click.command()
@click.option('--alignfile', help='input bam/sam file')
@click.option('--reffile', help='refrence fasta file (only for cram file)', default=None)
@click.option('--sitesfile', help='include sites (1-based) in this file (not only), four columns: contig pos ref alt', default=None)
@click.option('--minbasequlity', help='Minimum base quality (default=0). Bases below the minimum quality will not be counted.', default=0)
@click.option('--minmapqulity', help='Minimum mapping quality (default=0). Bases below the minimum quality will not be counted.', default=0)
@click.option('--persite', help='print persite result to this file, default is None', default=None)
@click.option('--outfile', help='output file')
def main(alignfile, reffile, sitesfile, minbasequlity, minmapqulity, persite, outfile):
    """
    输入的--sitesfile包含四列（contig pos(1-based) ref alt），没header，空格或tab分割
    """
    aln = load_alnfile(alignfile, reffile)
    sitedf = load_sites(sitesfile)
    result = []
    for contig, pos, ref, alt in sitedf.values:
        for record in aln.pileup(contig=contig, start=pos-1, end=pos, min_base_quality=0, min_mapping_quality=0):
            if record.reference_pos == (pos - 1):
                alleles = [x.upper() for x in record.get_query_sequences()]
                names = record.get_query_names()
                baseqs = record.get_query_qualities()
                mapqs = record.get_mapping_qualities()
                assert len(names) == len(alleles)
                assert len(names) == len(baseqs)
                assert len(names) == len(mapqs)
                for readname, allele, baseq, mapq in zip(names, alleles, baseqs, mapqs):
                    if (baseq < minbasequlity) or (mapq < minmapqulity):
                        print(f'low qulity: {readname}')
                        continue
                    if allele == ref:
                        refcount = 1
                        altcount = 0
                    elif allele == alt:
                        refcount = 0
                        altcount = 1
                    else:
                        refcount = 0
                        altcount = 0
                        print(contig, pos, ref, alt, readname, allele)
                        continue
                    result.append([readname, contig, pos, refcount, altcount, baseq, mapq])
    rdf = pd.DataFrame(result, columns=['name', 'contig', 'pos', 'refcount', 'altcount', 'baseQ', 'mapQ'])
    if persite:
        rdf.to_csv(persite, sep='\t', index=False)
    rdf.groupby('name').agg({'refcount': 'sum', 'altcount': 'sum'}).reset_index().to_csv(outfile, sep='\t', index=False)


if __name__ == '__main__':
    main()




