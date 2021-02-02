import re
import typer
import allel
import numpy as np
import pandas as pd
from scipy import stats


def check_mendelian(gts: list):
    """
    判断是否符合孟德尔遗传定律
    gts为一维列表，顺序为父本、母本、子代
    如['0/0', '1/1', '0/1']
    """
    if './.' in gts:
        return 'uninformative'
    F_alleles = set(re.split(r'[|,/]', gts[0]))
    M_alleles = set(re.split(r'[|,/]', gts[1]))
    C_alleles = set(re.split(r'[|,/]', gts[2]))
    combs = []
    for F in F_alleles:
        for M in M_alleles:
            combs.append(set([F, M]))
    if C_alleles in combs:
        return 'consistent'
    else:
        return 'inconsistent'


def classify_parents(gts: list):
    """
    返回子代中来自父母方的allele，无法判断时返回nan
    gts为一维列表，顺序为父本、母本、子代
    如['0/0', '1/1', '0/1']
    """
    if './.' in set(gts):
        return np.nan, np.nan
    F_alleles = set(re.split(r'[|,/]', gts[0]))
    M_alleles = set(re.split(r'[|,/]', gts[1]))
    F_sp_alleles = F_alleles - M_alleles
    M_sp_alleles = M_alleles - F_alleles
    if len(F_sp_alleles) > 0: # 二倍体
        F_alle = list(F_sp_alleles)[0]
        M_alle = list(M_alleles - F_sp_alleles)[0]
        return F_alle, M_alle
    elif len(M_sp_alleles) > 0:
        M_alle = list(M_sp_alleles)[0]
        F_alle = list(F_alleles - M_sp_alleles)[0]
        return F_alle, M_alle
    else:
        return np.nan, np.nan


def main(vcffile: str = typer.Argument(..., help="包含父母本DNA和RNA的总vcf文件"),
         sampleinfo: str = typer.Argument(..., help="样本信息文件，每行一个杂交RNA样本，tab分割的五列分别为 杂交RNA样本在VCF中的ID，杂交样本的组织类型，杂交个体的DNA在VCF中的ID，父本DNA在VCF中的ID，母本DNA在VCF中的ID"),
         outfile: str = typer.Argument(..., help="输出文件名(非压缩格式，.tsv后缀)")):
    cols = ['CHROM', 'POS', 'REF', 'ALTs', 'CID_rna', 'REF_count', 'ALTs_count', 'DP', 'MajorAllele_index', 'MajorAllele_parent', 'Major_Proportion', 'AllelicBias', 'BinomialTest', 'FisherExact', 'F_GT', 'M_GT', 'DNA_GT', 'RNA_GT', 'Mendelian_DNA', 'Mendelian_RNA']
    with open(outfile, 'w') as f:
        f.write('\t'.join(cols))
    with open(sampleinfo) as f:
        for line in f:
            if line[0] == '#':
                continue
            CID_rna, tissue, CID_dna, FID_dna, MID_dna = line.strip().split('\t')
            print(CID_rna, tissue, CID_dna, FID_dna, MID_dna)
            samples_query = [CID_rna, CID_dna, FID_dna, MID_dna]
            callset = allel.read_vcf(vcffile, samples=samples_query,
                                     fields=['samples', 'variants/CHROM', 'variants/POS', 'variants/REF', 'variants/ALT', 
                                             'calldata/GT', 'calldata/AD'])
            CID_rna_index = callset['samples'].tolist().index(CID_rna)
            try:
                CID_dna_index = callset['samples'].tolist().index(CID_dna)
            except ValueError:
                print(f'{CID_rna} 对应的DNA样本ID在VCF文件中不存在，孟德尔错误的检测将使用RNA的Genotype去做')
                CID_dna_index = None
            FID_dna_index = callset['samples'].tolist().index(FID_dna)
            MID_dna_index = callset['samples'].tolist().index(MID_dna)
            # 从VCF中提取原始信息
            df = pd.DataFrame({'CHROM': callset['variants/CHROM'],
                       'POS': callset['variants/POS'],
                       'REF': callset['variants/REF'],
                       'REF_count': callset['calldata/AD'][:, CID_rna_index, 0]})
            df['CID_rna'] = CID_rna
            df['ALTs'] = pd.DataFrame(callset['variants/ALT']).apply(lambda x: ','.join(x).strip(','), axis=1)

            # 各ALT的DP之和为ALTs_count
            alts_counts = callset['calldata/AD'][:, CID_rna_index, 1:]
            alts_counts[alts_counts == -1] = 0
            df['ALTs_count'] = alts_counts.sum(axis=1)

            # DP = REF + ALTs
            df['DP'] = df['REF_count'] + df['ALTs_count']

            # major allele相关计算
            df[['MajorAllele_count']] = callset['calldata/AD'][:, CID_rna_index, :].max(axis=1)
            df[['MinorAllele_count']] = df['DP'] - df['MajorAllele_count']
            df['MajorAllele_index'] = callset['calldata/AD'][:, CID_rna_index, :].argmax(axis=1)

            # DNA genotype信息
            gt_array = allel.GenotypeArray(callset['calldata/GT']).to_gt().astype(str)
            df['F_GT'] = gt_array[:, FID_dna_index]
            df['M_GT'] = gt_array[:, MID_dna_index]
            df['RNA_GT'] = gt_array[:, CID_rna_index]
            if CID_dna_index != None:
                df['DNA_GT'] = gt_array[:, CID_dna_index]
            else:
                df['DNA_GT'] = np.nan

            # 判断来源 （父方还是母方高表达）
            df['ParentsAllele'] = df[['F_GT', 'M_GT', 'DNA_GT']].apply(classify_parents, axis=1)
            df['F_allele_index'] = df['ParentsAllele'].str[0].astype(float)
            df['M_allele_index'] = df['ParentsAllele'].str[1].astype(float)
            df['MajorAllele_parent'] = np.nan
            df.loc[df['MajorAllele_index'] == df['F_allele_index'], 'MajorAllele_parent'] = 'F'
            df.loc[df['MajorAllele_index'] == df['M_allele_index'], 'MajorAllele_parent'] = 'M'

            # 孟德尔错误
            df['Mendelian_DNA'] = np.nan
            df['Mendelian_RNA'] = np.nan
            if CID_dna != None:
                df['Mendelian_DNA'] = df[['F_GT', 'M_GT', 'DNA_GT']].apply(check_mendelian, axis=1)
            else:
                df['Mendelian_RNA'] = df[['F_GT', 'M_GT', 'RNA_GT']].apply(check_mendelian, axis=1)

            # 一些差异表达指标
            df['Major_Proportion'] = df['MajorAllele_count'] / df['DP']
            df['AllelicBias'] = df['Major_Proportion'].values - 0.5 # major_pop总是大于0.5

            # 假设检验
            df['BinomialTest'] = df[['MajorAllele_count', 'DP']].apply(lambda ay: stats.binom_test(ay[0], ay[1], 0.5, 'two-sided'), axis=1)
            df['FisherExact'] = df[['MajorAllele_count', 'MinorAllele_count', 'DP']].apply(lambda ay: stats.fisher_exact([[ay[0], ay[1]], [ay[2]/2, ay[2]/2]])[1], axis=1)

            # 输出
            cols = ['CHROM', 'POS', 'REF', 'ALTs', 'CID_rna',
                    'REF_count', 'ALTs_count', 'DP',
                    'MajorAllele_index', 'MajorAllele_parent',
                    'Major_Proportion', 'AllelicBias', 'BinomialTest', 'FisherExact',
                    'F_GT', 'M_GT', 'DNA_GT', 'RNA_GT', 'Mendelian_DNA', 'Mendelian_RNA']
            df[cols].to_csv(outfile, sep='\t', index=False, mode='a', header=False)


if __name__ == '__main__':
    typer.run(main)