#!/usr/bin/env python3
from cyvcf2 import VCF
import argparse
import pandas as pd
import numpy as np

def main():
    parser = argparse.ArgumentParser(description="Count genotype categories using gt_types.")
    parser.add_argument("-i", "--input", required=True, help="Input VCF/VCF.gz")
    parser.add_argument("-o", "--output", required=True, help="Output TSV")
    args = parser.parse_args()

    vcf = VCF(args.input)
    samples = vcf.samples
    n_samples = len(samples)

    # 初始化计数数组（numpy 更快）
    het = np.zeros(n_samples, dtype=int)
    hom_ref = np.zeros(n_samples, dtype=int)
    hom_alt = np.zeros(n_samples, dtype=int)

    for variant in vcf:
        gt = variant.gt_types  # numpy array of 0,1,2,3

        # 直接按值分类计数
        hom_ref += (gt == 0)
        het     += (gt == 1)
        hom_alt += (gt == 3)
        # (gt == 2) 是 UNKNOWN，不计入

    # 计算 het ratio
    total = het + hom_ref + hom_alt
    het_ratio = np.divide(het, total, out=np.zeros_like(het, dtype=float), where=(total > 0))

    # 创建 DataFrame
    df = pd.DataFrame({
        "sample": samples,
        "het": het,
        "hom_ref": hom_ref,
        "hom_alt": hom_alt,
        "het_ratio": het_ratio
    })

    df.to_csv(args.output, sep="\t", index=False)
    print("Done. Results written to:", args.output)

if __name__ == "__main__":
    main()
