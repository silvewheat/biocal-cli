import sys
import pandas as pd

length_chromsome = 100000000
indfile = sys.argv[1]
df = pd.read_csv(indfile, sep='\t')
df['length'] = df['end'] - df['start'] + 1
gdf = df.groupby(['hapID', 'sourceID']).sum()['length'].reset_index()
gdf['percent'] = gdf['length'] / length_chromsome
gdf.to_csv('total_intro_per_ind.tsv', sep='\t', index=False)
