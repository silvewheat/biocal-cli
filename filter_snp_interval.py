# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 15:09:02 2020

@author: YudongCai

@Email: yudongcai216@gmail.com
"""

import click


@click.command()
@click.option('--invcf')
@click.option('--outvcf')
@click.option('--interval', type=int)
def main(invcf, outvcf, interval):
    print(f'select snp with a {interval} interval')
    with open(invcf) as f1, open(outvcf, 'w') as f2:
        for line in f1:
            if line[0] != '#':
                pos = int(line.split()[1])
                try:
                    dist = pos - lastpos
                except NameError:
                    lastpos = pos
                    f2.write(line)
                else:
                    lastpos = pos
                    if dist >= interval:
                        f2.write(line)
            else:
                f2.write(line)


if __name__ == '__main__':
    main()
