# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 15:24:16 2019

@author: YudongCai

@Email: yudongcai216@gmail.com
"""

import click
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC


@click.command()
@click.option('--infile', help='fasta file (mRNA)')
@click.option('--outfile', help='output file')
def main(infile, outfile):
    with open(outfile, 'w') as f:
        f.write(f'index\tlen_ORF\tlen_total\tmaxORF\n')
        for seq_record in SeqIO.parse(infile, 'fasta', IUPAC.unambiguous_dna):
            len_ORF = len(seq_record.translate(to_stop=True))
            len_total = len(seq_record)
            #Create three reading frames in forward direction, offset 0, 1, 2
            readingframes = [Seq.translate(seq_record.seq[i:], table='Standard', stop_symbol='*', to_stop=False, cds=False) for i in range(3)]
            results = []
            for frame in readingframes:
                for peptide in frame.split('*'): #Split translation over stopcodons
                    results.append(len(peptide))
            maxORF = max(results)

            f.write(f'{seq_record.id}\t{len_ORF}\t{len_total}\t{maxORF}\n')





if __name__ == '__main__':
    main()
