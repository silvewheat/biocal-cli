# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 15:20:26 2019

@author: YudongCai

@Email: yudongcai216@gmail.com
"""

import gzip
import click
import numpy as np
#import pandas as pd
import pyBigWig





@click.command()
@click.option('--bwfile', help='bigwig file for PhyloP')
@click.option('--beaglefile', help='genotype likelihoods (2) in beagle format (produced by angsd)')
@click.option('--majorfile', help='major allele in MAF, produced by maf_majorAllele.py')
@click.option('--outfile', help='output file name')
def main(bwfile, beaglefile, majorfile, outfile):
    """
    计算genetic load
    方法参考：
    Librado, P. et al. Ancient genomic changes associated with domestication of the horse. Science 356, 442–445 (2017).
    """
    print(__doc__)
    code2nuc = {'0': 'A',
                '1': 'C',
                '2': 'G',
                '3': 'T'}
#    nuc2code = {y: x for x, y in code2nuc.items()}
    majors = {}
    with gzip.open(majorfile, 'rb') as f:
        for line in f:
            tline = line.decode().strip().split()
            majors[tline[0]] = tline[1]
    markers = set(majors.keys())
    print('majorfile loaded')
    bw = pyBigWig.open(bwfile)
    with open(outfile, 'w') as f2, gzip.open(beaglefile, 'rb') as f1:
        header = f1.readline()
        for line in f1:
            tline = line.decode().strip().split()
            marker, allele1, allele2 = tline[:3]
#            print(f'in {marker}')
            if marker in markers:
                print(f'overlap {marker}')
                allele1 = code2nuc[allele1]
                chrom, loc = marker.split('_')
                phylop = bw.values(chrom, int(loc)-1, int(loc))[0]
                print(phylop)
                if phylop >= 1.5:
                    ref = majors[marker]
                    likelihoods = [float(i) if i != '0.333333' else np.nan for i in tline[3:]] # 0.333333为缺失
                    if ref == allele1:
                        gloads = phylop * np.array(likelihoods[2::3])
                    elif ref == allele2:
                        gloads = phylop * np.array(likelihoods[0::3])
                    else:
                        print(f'continue {marker}')
                        continue
                    print(f'oline {marker}')
                    olines = [chrom, loc] + [f'{x:.6f}' for x in gloads]
                    f2.write('\t'.join(olines) + '\n')


if __name__ == '__main__':
    main()






