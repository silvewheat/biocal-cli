# -*- coding: utf-8 -*-
"""
Created on Wed Oct 17 11:23:30 2018

@author: YudongCai

@Email: yudongcai216@gmail.com
"""


import click
import pysam
import numpy as np
import pandas as pd



@click.command()
@click.option('--vcffile', help='indexed vcf(.gz)/bcf file')
@click.option('--samplefile', help='query sample IDs, one ID per row')
@click.option('--outprefix', help='out file prefix')
def main(vcffile, samplefile, outprefix):
    """
    calssfy amd count the major allele in query samples
    """
    var = pysam.VariantFile(vcffile)
    samples = [x.strip() for x in open(samplefile).readlines()]
    var.subset_samples(samples)
    majors = []
    chrom = []
    pos = []
    counts = []
    for record in var:
        gts = np.array([x['GT'] for x in record.samples.values()])
        alleles = record.alleles
        count = []
        for n, allele in enumerate(alleles):
            count.append(np.sum(gts == n))
        count = np.array(count)
        majorallele = alleles[count.argmax()]
        majors.append(majorallele)
        chrom.append(record.chrom)
        pos.append(record.pos)
        counts.append(count.max())
    df = pd.DataFrame({'chrom': chrom, 'pos': pos, 'majorAllele': majors, 'alleleCount': counts})
    df[['chrom', 'pos', 'majorAllele', 'alleleCount']].to_csv(f'{outprefix}.tsv.gz', index=False, sep='\t', compression='gzip')


if __name__ == '__main__':
    main()
