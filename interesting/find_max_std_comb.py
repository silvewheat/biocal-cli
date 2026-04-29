import typer
import numpy as np
from itertools import combinations


def main(infile: str = typer.Option(..., help="Input file, one column"),
         valuecol: int = typer.Option(..., help="第几列是值(1-base)"),
         combsize: int = typer.Option(..., help="取多少行"),
         outfile: str = typer.Option(..., help="输出文件")):
    valuecol -= 1
    values = []
    with open(infile) as f:
        for line in f:
            items = line.strip().split()
            values.append(float(items[valuecol]))
    indexs = np.argsort(values).tolist()
    values = sorted(values)
    nlines = len(values)
    max_std = 0
    for min_n in range(1, combsize):
        comb = values[ :min_n] + values[-(combsize-min_n): ]
        comb_index = indexs[ :min_n] + indexs[-(combsize-min_n): ]
        std = np.std(comb)
        if std > max_std:
            max_std = std
            best_comb = comb
            best_comb_index = comb_index

    best_comb_index = set(best_comb_index)
    with open(infile) as f1, open(outfile, 'w') as f2:
        for nline, line in enumerate(f1):
            if nline in best_comb_index:
                f2.write(line)

if __name__ == "__main__":
    typer.run(main)
