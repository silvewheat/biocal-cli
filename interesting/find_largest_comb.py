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
    values = np.array(values)
    nlines = len(values)
    comb_best = ()
    std_largest = 0
    for comb_index in combinations(np.arange(nlines), combsize):
        std = np.std(values[list(comb_index)])
        if std > std_largest:
            std_largest = std
            comb_best = comb_index
    comb_best = set(comb_best)
    with open(infile) as f1, open(outfile, 'w') as f2:
        for nline, line in enumerate(f1):
            if nline in comb_best:
                f2.write(line)

if __name__ == "__main__":
    typer.run(main)
