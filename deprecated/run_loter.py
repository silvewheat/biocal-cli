# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 10:26:47 2020

@author: YudongCai

@Email: yudongcai216@gmail.com
"""

import pickle
import sys
import click
import allel
import numpy as np
import pandas as pd
import loter.locanc.local_ancestry as lc


def vcf2npy(vcffile, samples):
    callset = allel.read_vcf(vcffile, samples=samples,
                             fields=['samples', 'variants/CHROM', 'variants/POS', 'calldata/GT'])
    haplotypes_1 = callset['calldata/GT'][:,:,0]
    haplotypes_2 = callset['calldata/GT'][:,:,1]

    m, n = haplotypes_1.shape
    mat_haplo = np.empty((2*n, m))
    mat_haplo[::2] = haplotypes_1.T
    mat_haplo[1::2] = haplotypes_2.T

    keep_samples = callset['samples']
    return mat_haplo.astype(np.uint8), keep_samples, callset['variants/CHROM'], callset['variants/POS']


@click.command()
@click.option('--vcffile')
@click.option('--refpops', multiple=True, help='多个ref申明多次')
@click.option('--querypop')
@click.option('--thread', default=1, type=int)
@click.option('--outprefix')
@click.option('--plot-heatmap', is_flag=True, default=False, help='plot heatmap')
def main(vcffile, refpops, querypop, thread, outprefix, plot_heatmap):
    """
    run loter
    querypop is mixed between refpop1, refpop2
    """
    # load query
    samples = [x.strip() for x in open(querypop)]
    H_query, ID_query, chrlist, poslist = vcf2npy(vcffile, samples)
    print('query loaded.')
    print(f'{len(ID_query)} individuals, {H_query.shape[1]} sites.')
    chrom = chrlist[0]

    # load ref
    H_refs = []
    for refpop in refpops:
        samples = [x.strip() for x in open(refpop)]
        H_ref, IDs, chrs, sites = vcf2npy(vcffile, samples)
        print(f'ref {refpop} loaded.')
        print(f'{len(IDs)} individuals, {H_ref.shape[1]} sites.')
        H_refs.append(H_ref)
    del(H_ref)
    del(sites)
    del(chrs)
    del(chrlist)

    # calculate
    print('runing loter')
    if len(refpops) == 2:
        res_loter = lc.loter_smooth(l_H=H_refs, h_adm=H_query,
                                    range_lambda=np.arange(1.5, 5.5, 0.5),
                                    threshold=0.90, rate_vote=0.5, nb_bagging=20,
                                    num_threads=thread)
    elif len(refpops) > 2:
        print('多于两个refpop还没写好')
        res_impute, res_no_impute = lc.loter_local_ancestry(l_H=H_refs, h_adm=H_query, num_threads=thread)
        print('我先把结果存起来')
        with open(f'{outprefix}_impute.pkl', 'wb') as f:
            pickle.dump(res_impute, f)
        with open(f'{outprefix}_no_impute.pkl', 'wb') as f:
            pickle.dump(res_no_impute, f)
        print('存完了，我先退出了')
        sys.exit()
        
    print('loter finished')
    del(H_query)
    np.save(f'{outprefix}', res_loter)
    hapIDs = [[f'{i}_1', f'{i}_2'] for i in ID_query]
    hapIDs = [i for j in hapIDs for i in j]
    df = pd.DataFrame(res_loter.T, columns=hapIDs)
    df['chr'] = chrom
    df['pos'] = poslist
    cols = ['chr', 'pos'] + hapIDs
    # write the per SNP classification
    print('write perSNP output')
    df[cols].to_csv(f'{outprefix}_perSNPs.tsv.gz', sep='\t', index=False, compression='gzip')

    # plot heatmap
    if plot_heatmap:
        import matplotlib.pyplot as plt
        import sebaorn as sns
        print('plotting heatmap')
        f, ax = plt.subplots(figsize=(20, 10))
        df.index = df['pos']
        ax = sns.heatmap(df.T, yticklabels=1)
        plt.savefig(f'{outprefix}.jpg', dpi=300)
        plt.close()

    # merge the consecutive rows with same values
    print('merge introgressed segments in pop level')
    cumindex = (df[hapIDs] != df[hapIDs].shift()).apply(max, axis=1).cumsum()
    fundict = {'pos': [min, max]}
    fundict.update(dict(zip(hapIDs, [max]*len(hapIDs)))) # same for 'max' and 'min'
    mdf1 = df.groupby(cumindex).agg(fundict)
    mdf1.columns = ['start', 'end'] + hapIDs
    mdf1['chr'] = chrom
    cols = ['chr', 'start', 'end'] + hapIDs
    mdf1[cols].to_csv(f'{outprefix}_popSeg.tsv.gz', sep='\t', index=False, compression='gzip')
    del(mdf1)

    # merge the consecutive rows for each indiv
    print('merge introgressed segments in indiv level')
    df = df.melt(id_vars=['chr', 'pos'], value_vars=hapIDs,
                 var_name='hapID', value_name='sourceID')
    cumindex = (df[['hapID', 'sourceID']] != df[['hapID', 'sourceID']].shift()).apply(max, axis=1).cumsum()
    fundict = {'chr': max, 'pos': [min, max], 'hapID': max, 'sourceID': max}
    mdf2 = df.groupby(cumindex).agg(fundict)
    mdf2.columns = ['chr', 'start', 'end', 'hapID', 'sourceID']
    mdf2.to_csv(f'{outprefix}_indSeg.tsv.gz', sep='\t', index=False, compression='gzip')
    print('Done!')


if __name__ == '__main__':
    main()

