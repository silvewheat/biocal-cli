# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 14:25:51 2019

@author: YudongCai

@Email: yudongcai216@gmail.com
"""

import click
import numpy as np
import pysam



def count_allenum(gts):
    num_alleles = 0
    for gt in gts.flatten():
        if gt != None:
            num_alleles += 1
    return num_alleles

def count_refnum(gts, ref):
    return np.sum(gts.flatten() == ref)



@click.command()
@click.option('--vcffile', help='.vcf/.vcf.gz/.bcf')
@click.option('--samplefile', help='samples list')
@click.option('--outfile', help='outfile name')
def main(vcffile, samplefile, outfile):
    samples = [x.strip() for x in open(samplefile)]
    var = pysam.VariantFile(vcffile)
    var.subset_samples(samples)
    with open(outfile, 'w') as f:
        f.write('chrom\tpos\tREF_freq\trefCount\talleleCount\n')
        for record in var:
            ref = record.ref
            gts = np.array([x.alleles for x in record.samples.values()])
            refcount = count_refnum(gts, ref)
            allelecount = count_allenum(gts)
            if allelecount > 0:
                reffreq = refcount / allelecount
                f.write(f'{record.contig}\t{record.pos}\t{reffreq:.2f}\t{refcount}\t{allelecount}\n')
            else:
                f.write(f'{record.contig}\t{record.pos}\tnan\t{refcount}\t{allelecount}\n')

if __name__ == '__main__':
    main()


