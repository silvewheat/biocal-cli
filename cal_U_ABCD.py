# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 17:39:45 2020

@author: YudongCai

@Email: yudongcai216@gmail.com
"""

import click
import allel
import numpy as np


@click.command()
@click.option('--vcffile')
@click.option('--popa')
@click.option('--popb')
@click.option('--popc')
@click.option('--popd')
@click.option('--freqw', type=float)
@click.option('--freqx', type=float)
@click.option('--freqy', type=float)
@click.option('--freqz', type=float)
@click.option('--outfile')
@click.option('--outfreqfile', default=None)
def main(vcffile, popa, popb, popc, popd, freqw, freqx, freqy, freqz, outfile, outfreqfile):
    """
    U_A,B,C,D(w,x,y,z)
    A是非渗入群体，B是被渗入群体，C是渗入来源群体1, D是渗入来源群体2
    在窗口内A中频率小于w，B中大于x，C中大于y，D中小于z的SNP位点数即为U_A,B,C,D(w,x,y,z)
    详见：Signatures of Archaic Adaptive Introgression in Present-Day Human Populations
    """
    popA = [x.strip() for x in open(popa)]
    popB = [x.strip() for x in open(popb)]
    popC = [x.strip() for x in open(popc)]
    popD = [x.strip() for x in open(popd)]
    callset_C = allel.read_vcf(vcffile, samples=popC,
                             fields=['samples', 'variants/CHROM', 'variants/POS', 'calldata/GT'])
    gt_C = allel.GenotypeArray(callset_C['calldata/GT'])
    ac_C = gt_C.count_alleles()
    af_C = ac_C.to_frequencies()
    site_selection = np.sum(af_C >= freqy, axis=1) > 0 # 只保留C群体中频率大于y的位点
    pos = callset_C['variants/POS'][site_selection]
    chroms = callset_C['variants/CHROM'][site_selection]
    allel_selection = af_C >= freqy # 筛选allele的编号，以C群体中频率最大的allel为准（包含了site_selection的内容）
    af_C = af_C[allel_selection]
    del(callset_C)
    del(gt_C)
    del(ac_C)

    callset_A = allel.read_vcf(vcffile, samples=popA,
                             fields=['calldata/GT'])
    af_A = allel.GenotypeArray(callset_A['calldata/GT']).count_alleles().to_frequencies()[allel_selection]
    del(callset_A)

    callset_B = allel.read_vcf(vcffile, samples=popB,
                             fields=['calldata/GT'])
    af_B = allel.GenotypeArray(callset_B['calldata/GT']).count_alleles().to_frequencies()[allel_selection]
    del(callset_B)

    callset_D = allel.read_vcf(vcffile, samples=popD,
                             fields=['calldata/GT'])
    af_D = allel.GenotypeArray(callset_D['calldata/GT']).count_alleles().to_frequencies()[allel_selection]
    del(callset_D)

    Usites_selection = (af_A <= freqw) & (af_B >= freqx) & (af_C >= freqy) & (af_D <= freqz)
    U_chroms = chroms[Usites_selection]
    U_pos = pos[Usites_selection]

    with open(outfile, 'w') as f:
        for chrom, pos in zip(U_chroms, U_pos):
            f.write(f'{chrom}\t{pos}\n')
    with open(outfreqfile, 'w') as f:
        f.write('chrom\tpos\tfreqA\tfreqB\tfreqC\tfreqD\n')
        for chrom, pos, freqA, freqB, freqC, freqD in zip(chroms, pos, af_A, af_B, af_C, af_D):
            f.write(f'{chrom}\t{pos}\t{freqA:.3f}\t{freqB:.3f}\t{freqC:.3f}\t{freqD:.3f}\n')

if __name__ == '__main__':
    main()



