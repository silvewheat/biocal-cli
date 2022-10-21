# -*- coding: utf-8 -*-
"""
Created on Mon Dec  3 10:46:11 2018

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

def load_region(regionfile):
    regions = []
    with open(regionfile) as f:
        for line in f:
            if line[0] != '#':
                sline = line.strip().split()
                regions.append([sline[0], int(sline[1]), int(sline[2])])
    return regions


def load_sites(sitesfile):
    regions = []
    with open(sitesfile) as f:
        for line in f:
            if line[0] != '#':
                sline = line.strip().split()
                regions.append([sline[0], int(sline[1])-1, int(sline[1])]) # 1-base to 0-base
    return regions





@click.command()
@click.option('--alignfile', help='input bam/sam file')
@click.option('--reffile', help='refrence fasta file (only for cram file)', default=None)
@click.option('--regionfile', help='restrict to regions in this file. three columns: chrom start end', default=None)
@click.option('--sitesfile', help='include sites (1-based) in this file (not only), two columns: chrom pos', default=None)
@click.option('--minbasequlity', help='Minimum base quality (default=0). Bases below the minimum quality will not be counted.', default=0)
@click.option('--outfile', help='output file')
def main(alignfile, reffile, regionfile, sitesfile, minbasequlity, outfile):
    aln = load_alnfile(alignfile, reffile)
    nucs = {'A', 'C', 'G', 'T', 'a', 'c', 'g', 't'}
    results = []
    if regionfile:
        regions = load_region(regionfile)
        for region in regions:
            print(region)
            chrom = region[0]
            for record in aln.pileup(*region, min_base_quality=minbasequlity):
                try:
                    if record.get_num_aligned() > 0:
                        alleles = [x.upper() for x in record.get_query_sequences() if x in nucs]
                        if len(np.unique(alleles)) > 1:
                            pos = record.reference_pos
                            depth = len(alleles)
                            for allele, count in Counter(alleles).items():
                                results.append([chrom, pos+1, allele, count, depth])
                except AssertionError:
                    print(record.pos, record.get_num_aligned())
                    pass
    if sitesfile:
        regions = load_sites(sitesfile)
        for region in regions:
            chrom = region[0]
            loc = region[1]
            for record in aln.pileup(*region, min_base_quality=minbasequlity):
                pos = record.reference_pos # 0-based
                if pos == loc:
                    if record.get_num_aligned() > 0:
                        alleles = [x.upper() for x in record.get_query_sequences() if x in nucs]
                        depth = len(alleles)
                        for allele, count in Counter(alleles).items():
                            results.append([chrom, pos+1, allele, count, depth])
    results = pd.DataFrame(results, columns=['chrom', 'pos', 'allele', 'count', 'depth'])
    comp = 'gzip' if outfile.split('.')[-1] == 'gz' else None
    results.drop
    results.sort_values(['chrom', 'pos'], inplace=True)
#    results.to_csv(outfile, sep='\t', index=False, compression=comp)
    results = results.pivot_table(index=['chrom', 'pos'], columns=['allele'], values=['count'])
    results.columns = results.columns.get_level_values(1)
    results.reset_index(inplace=True)
    results.to_csv(outfile, sep='\t', compression=comp, index=False)


if __name__ == '__main__':
    main()






