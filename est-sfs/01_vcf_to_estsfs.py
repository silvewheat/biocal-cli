# -*- coding: utf-8 -*-
"""
Created on 2022 10-14

@author: Yudongcai

@Email: yudong_cai@163.com
"""
import re
import typer
import numpy as np
from cyvcf2 import VCF
from collections import Counter, defaultdict


def convert_gts(gt_bases):
    gt_split = re.compile(r'[/|]')
    bases = []
    for base in gt_bases:
        bases.extend(gt_split.split(base))
    return bases


def main(vcffile: str = typer.Argument(..., help="input vcf file"),
         focalsamples: str = typer.Argument(..., help="sample list for focal samples"),
         outgroup1: str = typer.Argument(..., help="sample list for outgroup1"),
         outgroup2: str = typer.Argument(..., help="sample list for outgroup2"),
         outgroup3: str = typer.Argument(..., help="sample list for outgroup3"),
         outprefix: str = typer.Argument(..., help="output prefix")):
    focal_samples = [x.strip() for x in open(focalsamples)]
    outgroup1_samples = [x.strip() for x in open(outgroup1)]
    outgroup2_samples = [x.strip() for x in open(outgroup2)]
    outgroup3_samples = [x.strip() for x in open(outgroup3)]
    samples = focal_samples + outgroup1_samples + outgroup2_samples + outgroup3_samples
    print(f'focal samples: {len(focal_samples)}\noutgroup1: {len(outgroup1_samples)}\noutgroup2: {len(outgroup2_samples)}\noutgroup3: {len(outgroup3_samples)}')
    
    with open(f'{outprefix}_siteInfo.tsv', 'w') as f1, open(f'{outprefix}_datafile', 'w') as f2:
        base2index = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        f1.write('CHROM\tPOS\tREF\tALT\tmajorAllele\tminorAllele\n')
        vcf = VCF(vcffile, gts012=True, samples=samples)
        focal_selection = [True if x in focal_samples else False for x in vcf.samples]
        outgroup1_selection = [True if x in outgroup1_samples else False for x in vcf.samples]
        outgroup2_selection = [True if x in outgroup2_samples else False for x in vcf.samples]
        outgroup3_selection = [True if x in outgroup3_samples else False for x in vcf.samples]
        outgroup_selections = (outgroup1_selection, outgroup2_selection, outgroup3_selection)
        for variant in vcf:
            alleles = [variant.REF] + variant.ALT
            f1.write(f'{variant.CHROM}\t{variant.POS}\t{variant.REF}\t' + ','.join(variant.ALT) + '\t')
            counter_gts_focal = Counter(convert_gts(variant.gt_bases[focal_selection]))
            major_allele = counter_gts_focal.most_common()[0][0]
            try:
                minor_allele = counter_gts_focal.most_common()[1][0]
            except IndexError:
                minor_allele = list(set(alleles) - set(major_allele))[0]
            f1.write(f'{major_allele}\t{minor_allele}\n')
            f2.write(f"{counter_gts_focal.get('A', 0)},{counter_gts_focal.get('C', 0)},{counter_gts_focal.get('G', 0)},{counter_gts_focal.get('T', 0)}")

            for selection in outgroup_selections:
                counts = ['0', '0', '0', '0'] # A C G T
                counter_gts = Counter(convert_gts(variant.gt_bases[selection])).most_common()
                first_base, first_count = counter_gts[0]
                try:
                    second_base, second_count = counter_gts[1]
                except IndexError:
                    second_count = 0
                # 两种allele数量相等时按缺失处理
                if (first_count > second_count) and (first_base != '.'):
                    counts[base2index[first_base]] = '1'
                f2.write('\t'+','.join(counts))
            f2.write('\n')

if __name__ == '__main__':
    typer.run(main)