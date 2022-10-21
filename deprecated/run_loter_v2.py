# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 11:50:05 2020

@author: YudongCai

@Email: yudongcai216@gmail.com
"""

import click
import allel
import numpy as np
import pandas as pd
import loter.locanc.local_ancestry as lc


@click.command()
@click.option('--vcffile', help='只包含一条染色体的vcf文件')
@click.option('--chrom', help='输入文件的染色体ID，只用于输出，不会对vcf进行索引')
@click.option('--refpop1')
@click.option('--refpop2')
@click.option('--querypop')
@click.option('--thread', default=1, type=int)
@click.option('--outprefix')
@click.option('--plot-heatmap', is_flag=True, default=False, help='画热图')
def main(vcffile, chrom, refpop1, refpop2, querypop, thread, outprefix, plot_heatmap):
    """
    run loter
    querypop is mixed between refpop1, refpop2
    """
    # load vcf
    qsamples = [x.strip() for x in open(querypop)]
    r1samples = [x.strip() for x in open(refpop1)]
    r2samples = [x.strip() for x in open(refpop2)]
    allsamples = list(set(qsamples + r1samples + r2samples))
    callset = allel.read_vcf(vcffile, samples=allsamples,
                             fields=['samples', 'variants/POS', 'calldata/GT'])
    print('vcf loaded')
    poslist = callset['variants/POS']
    haplotypes_1 = callset['calldata/GT'][:,:,0]
    haplotypes_2 = callset['calldata/GT'][:,:,1]

    callset_samples = np.array(callset['samples'])

    del(callset)

    q_callset_index = [int(np.where(callset_samples == s)[0]) for s in qsamples]
    r1_callset_index = [int(np.where(callset_samples == s)[0]) for s in r1samples]
    r2_callset_index = [int(np.where(callset_samples == s)[0]) for s in r2samples]

    q_out_samples = callset_samples[q_callset_index]
    r1_out_samples = callset_samples[r1_callset_index]
    r2_out_samples = callset_samples[r2_callset_index]

    nsite, _ = haplotypes_1.shape

    print(f'nsites: {nsite}')
    print(f'n_query: {len(q_out_samples)}, n_ref1: {len(r1_out_samples)}, n_ref2: {len(r2_out_samples)}')

    H_query = np.empty((2*len(q_callset_index), nsite), dtype=np.uint8)
    H_query[::2] = haplotypes_1[:, q_callset_index].T
    H_query[1::2] = haplotypes_2[:, q_callset_index].T

    H_ref1 = np.empty((2*len(r1_callset_index), nsite), dtype=np.uint8)
    H_ref1[::2] = haplotypes_1[:, r1_callset_index].T
    H_ref1[1::2] = haplotypes_2[:, r1_callset_index].T

    H_ref2 = np.empty((2*len(r2_callset_index), nsite), dtype=np.uint8)
    H_ref2[::2] = haplotypes_1[:, r2_callset_index].T
    H_ref2[1::2] = haplotypes_2[:, r2_callset_index].T

    del(haplotypes_1)
    del(haplotypes_2)
    print('numpy array created')

    # calculate
    print('runing loter')
    res_loter = lc.loter_smooth(l_H=[H_ref1, H_ref2], h_adm=H_query,
                                range_lambda=np.arange(1.5, 5.5, 0.5),
                                threshold=0.90, rate_vote=0.5, nb_bagging=20,
                                num_threads=thread)
    print('loter finished')
    del(H_query)
    del(H_ref1)
    del(H_ref2)
    np.save(f'{outprefix}', res_loter)
    hapIDs = [[f'{i}_1', f'{i}_2'] for i in q_out_samples]
    hapIDs = [i for j in hapIDs for i in j]
    df = pd.DataFrame(res_loter.T, columns=hapIDs)
    df['chr'] = chrom
    df['pos'] = poslist
    cols = ['chr', 'pos'] + hapIDs
    # write the per SNP classification
    print('write perSNP output')
    df[cols].to_csv(f'{outprefix}_perSNPs.tsv.gz', sep='\t', index=False)
    # plot heatmap
    if plot_heatmap:
        import matplotlib.pyplot as plt
        import seaborn as sns
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
    fundict.update(dict(zip(hapIDs, [max]*len(hapIDs)))) # 这儿用min或者max都一样
    mdf1 = df.groupby(cumindex).agg(fundict)
    mdf1.columns = ['start', 'end'] + hapIDs
    mdf1['chr'] = chrom
    cols = ['chr', 'start', 'end'] + hapIDs
    mdf1[cols].to_csv(f'{outprefix}_popSeg.tsv.gz', sep='\t', index=False)
    del(mdf1)
    # merge the consecutive rows for each indiv
    print('merge introgressed segments in indiv level')
    df = df.melt(id_vars=['chr', 'pos'], value_vars=hapIDs,
                 var_name='hapID', value_name='sourceID')
    cumindex = (df[['hapID', 'sourceID']] != df[['hapID', 'sourceID']].shift()).apply(max, axis=1).cumsum()
    fundict = {'chr': max, 'pos': [min, max], 'hapID': max, 'sourceID': max}
    mdf2 = df.groupby(cumindex).agg(fundict)
    mdf2.columns = ['chr', 'start', 'end', 'hapID', 'sourceID']
    mdf2.to_csv(f'{outprefix}_indSeg.tsv.gz', sep='\t', index=False)
    print('Done!')


if __name__ == '__main__':
    main()
