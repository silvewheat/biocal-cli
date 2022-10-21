# -*- coding: utf-8 -*-
"""
Created on Sun Jul 21 14:53:57 2019

@author: YudongCai

@Email: yudongcai216@gmail.com
"""

import os
os.environ['KERAS_BACKEND']='tensorflow'
import tensorflow as tf
import click
import khmer
import numpy as np
import pandas as pd
from itertools import product
from keras.models import load_model


def tf_no_warning():
    """
    Make Tensorflow less verbose
    """
    try:
        tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)
        os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

    except ImportError:
        pass

def kmerproduce(infile, k, outfile):
    kmerlist = []
    for kmer in product('ACGT', repeat=k):
        kmerlist.append(''.join(kmer))

    df = []
    seqids = []
    for seq in khmer.ReadParser(infile):
        seqids.append(seq.name)
        df.append([])
        counts = khmer.Counttable(k, 1e6, 4)
        counts.set_use_bigcount(True)
        counts.consume(seq.sequence.upper())
        for kmer in kmerlist:
            df[-1].append(counts.get(kmer))

    df = pd.DataFrame(df, columns=kmerlist, index=seqids)
    df.T.reset_index().to_csv(outfile, index=False, sep='\t')


def load_features(outdir, outprefix):
    alldata = []
    for k in range(1, 6):
        print(f'k = {k}')
        file = os.path.join(outdir, f'{outprefix}_{k}mer.tsv.gz')
        df = pd.read_csv(file, sep='\t').set_index('index')
        df = df / df.sum()
        data = df.T.values
        alldata.append(data)
    return np.hstack(alldata)


def load_seqid(infile):
    seqids = []
    with open(infile) as f:
        for line in f:
            if line[0] == '>':
                seqids.append(line.strip().split()[0][1:])
    return seqids



@click.command()
@click.option('--infile', help='输入的fasta文件路径')
@click.option('--modelfile', help='模型路径')
@click.option('--outdir', help='输出文件目录')
@click.option('--outprefix', help='输出文件名前缀')
def main(infile, modelfile, outdir, outprefix):
    tf_no_warning()
    print('kmer producing...')
    for k in range(1, 6):
        print(f'k = {k}')
        outfile = os.path.join(outdir, f'{outprefix}_{k}mer.tsv.gz')
        kmerproduce(infile, k, outfile)

    print('loading features...')
    data = load_features(outdir, outprefix)

    print('loading model')
    model = load_model(modelfile)

    print('predicting')
    labels = model.predict_classes(data).flatten()
    seqids = load_seqid(infile)
    df = pd.DataFrame({'seqid': seqids, 'label': labels})
    outfile = os.path.join(outdir, f'{outprefix}_result.tsv.gz')
    df.to_csv(outfile, sep='\t', index=False)
    print('Completed!')


if __name__ == '__main__':
    main()
