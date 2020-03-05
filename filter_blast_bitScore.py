# -*- coding: utf-8 -*-
"""
Created on Wed Dec 25 15:07:13 2019

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
    对于每个query，选取bit score最大的相应结果，统计结果中会输出qcovhsp
    blast输入格式要求outfmt 7 std，外加一个qcovhsp就可以，可以再增加其他的字段
    比如为 -outfmt '7 std qcovs qcovhsp qcovus'
    """
    records = []
    scores = []
    with open(blastout) as f1, open(f'{outprefix}.bitScore.tsv', 'w') as f2, open(f'{outprefix}.stat.tsv', 'w') as f3:
        f3.write('query\tident\tqcovhsp\taln_length\n')
        for line in f1:
            if line[0] != '#':
                records.append(line)
                score_index = header.index('bit score')
                score = float(line.strip().split()[score_index])
                scores.append(score)
            else:
                f2.write(line)
                if line[2:7] == 'Query':
#                    print(line.strip('#'), end='')
                    query = line.strip().split()[2]
                if line[2:8] == 'Fields':
                    header = [x.strip() for x in line[10:].split(',')]
 #                   print(header)
                if records:
                    scores = np.array(scores)
                    max_score = scores.max()
                    for record, ismax in zip(records, scores==max_score):
                        if ismax:
                            f2.write(record)
                            tline = record.split()
                            ident = float(tline[header.index('% identity')])
                            length = int(tline[header.index('alignment length')])
                            qcovhsp = int(tline[header.index('% query coverage per hsp')])
                            f3.write(f'{query}\t{ident}\t{qcovhsp}\t{length}\n')
                    records = []
                    scores = []


if __name__ == '__main__':
    main()

