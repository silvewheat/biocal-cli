# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 11:20:06 2019

@author: YudongCai

@Email: yudongcai216@gmail.com
"""

import click
import numpy as np


@click.command()
@click.option('--blastout', help='blast输出结果, 最后一行#开头')
@click.option('--outprefix', help='过滤后的blast结果')
def main(blastout, outprefix):
    """
    对于每个query，选取query coverage per subject(qcovs)最大的相应结果
    blast输入格式要求outfmt 7 std，外加一个qcovs就可以，可以再增加其他的字段
    比如为 -outfmt '7 std qcovs qcovhsp qcovus'
    """
    records = []
    qcovs = []
    with open(blastout) as f1, open(f'{outprefix}.qcovs.tsv', 'w') as f2, open(f'{outprefix}.qstat.tsv', 'w') as f3:
        for line in f1:
            if line[0] != '#':
                records.append(line)
                qcov_index = header.index('% query coverage per subject')
                qcov = int(line.strip().split()[qcov_index])
                qcovs.append(qcov)
            else:
                f2.write(line)
                if line[2:7] == 'Query':
                    print(line.strip('#'), end='')
                    query = line.strip().split()[-1]
                if line[2:8] == 'Fields':
                    header = [x.strip() for x in line[10:].split(',')]
                    print(header)
                if records:
                    qcovs = np.array(qcovs)
                    max_qcov = qcovs.max()
                    lens = []
                    idents = []
                    for record, ismax in zip(records, qcovs==max_qcov):
                        if ismax:
                            f2.write(record)
                            tline = record.split()
                            lens.append(int(tline[header.index('alignment length')]))
                            idents.append(float(tline[header.index('% identity')]))
                    ident = np.average(idents, weights=lens)
                    f3.write(f'{query}\t{ident}\t{max_qcov}\n')
                    records = []
                    qcovs = []


if __name__ == '__main__':
    main()

