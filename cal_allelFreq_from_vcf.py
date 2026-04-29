import typer
from cyvcf2 import VCF
import numpy as np
import gzip


def main(vcffile: str = typer.Option(..., help="input vcf file"),
         samples: str = typer.Option(..., help="query sample file, one sample ID per line"),
         outprefix: str = typer.Option(..., help="output prefix")):
    """
    stat allele freq from vcf file
    output file: '{outprefix}_allelFreq.tsv.gz'
    """
    query_samples = [x.strip() for x in open(samples)]
    with gzip.open(f'{outprefix}_allelFreq.tsv.gz', 'wb') as f:
        f.write('chrom\tpos\tref\talt\tn_homRef\tn_het\tn_homAlt\talt_allel_freq\tn_with_ref_support\tn_with_alt_support\tn_het_with_MoreThanOne_support\n'.encode())
        for variant in VCF(vcffile, gts012=True, samples=query_samples):
            n_homRef = np.sum(variant.gt_types == 0)
            n_het = np.sum(variant.gt_types == 1)
            n_homAlt = np.sum(variant.gt_types == 2)
            alt_allel_freq = (2*n_homAlt + n_het) / ((n_homRef+n_het+n_homAlt) * 2)
            n_with_ref_support = np.sum(variant.gt_ref_depths > 0)
            n_with_alt_support = np.sum(variant.gt_alt_depths > 0)
            n_het_with_MoreThanOne_support = np.sum((variant.gt_ref_depths > 0) & (variant.gt_alt_depths > 0))
            outline = '\t'.join([str(x) for x in [variant.CHROM, variant.POS, variant.REF, variant.ALT[0], n_homRef, n_het, n_homAlt, alt_allel_freq, n_with_ref_support, n_with_alt_support, n_het_with_MoreThanOne_support, '\n']])
            f.write(outline.encode())

if __name__ == '__main__':
    typer.run(main)