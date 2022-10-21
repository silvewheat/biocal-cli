# -*- coding: utf-8 -*-
"""
Created on Fri Sep 14 20:23:43 2018

@author: YudongCai

@Email: yudongcai216@gmail.com
"""


import click



@click.command()
@click.argument('fafile')
def main(fafile):
    """

    """
    with open(fafile) as f:
        length = 0
        lastname = None
        for line in f:
            if line[0] != '>':
                length += len(line.strip())
            else:
                name = line[1:].split()[0]
                if lastname:
                    print(f'{lastname}\t{length}')
                lastname = name
                length = 0
        print(f'{lastname}\t{length}\n')


if __name__ == '__main__':
    main()
