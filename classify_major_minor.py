# -*- coding: utf-8 -*-
"""
Created on Sat Dec  1 17:01:50 2018

@author: YudongCai

@Email: yudongcai216@gmail.com
"""

import os
import gzip
import click
from collections import Counter



@click.command()
@click.option('--varfile', help='input vcf(.gz)/bcf file')
@click.option('--querysamples', help='query sample list')
@click.option('--outfile', help='outfile is gzipped')
def main(varfile, querysamples, outfile):
    cmd = f"bcftools query -S {querysamples} -f '%CHROM\t%POS\t%REF\t%ALT\t[%TGT\t]\n' {varfile}"
    with gzip.open(outfile, 'wb') as f:
        header = 'chrom\tpos\tmajor\tminor\n'.encode()
        f.write(header)
        for line in os.popen(cmd):
            tline = line.strip().split()
            chrom, pos, ref, alt = tline[:4]
            gts = [x[0] for x in tline[4:]] + [x[-1] for x in tline[4:]]
            count = Counter(gts)
            if count.most_common()[0][0] == ref:
                f.write(f'{chrom}\t{pos}\t{ref}\t{alt}\n'.encode())
            elif count.most_common()[0][0] == alt:
                f.write(f'{chrom}\t{pos}\t{alt}\t{ref}\n'.encode())
            else:
                print(f'{chrom}:{pos} seems not biallelic')

if __name__ == '__main__':
    main()

