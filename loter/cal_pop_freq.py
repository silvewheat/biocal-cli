import sys
import pandas as pd
from collections import defaultdict
from collections import Counter

source_index = int(sys.argv[3])

infile = sys.argv[1]
df = pd.read_csv(infile, sep='\t')

groupfile = sys.argv[2]
group2inds = defaultdict(list)
with open(groupfile) as f: 
    for line in f: 
        tline = line.strip().split() 
        group2inds[tline[1]].append(tline[0])

odf = df[['chr', 'start', 'end']].copy()
for group, inds in group2inds.items(): 
    odf[group] = df[inds].apply(lambda x: Counter(x).get(source_index), axis=1) / len(inds)

odf.to_csv('total_intro_freq.tsv', sep='\t', index=False)
