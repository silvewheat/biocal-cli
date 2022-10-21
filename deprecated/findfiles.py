import os
import glob

with open('filelist', 'w') as f:
    for file in glob.glob('regions_*.tsv.gz'):
        if len(file.split('.')) == 3:
            outfile = file[:-7] + '.mdf.tsv.gz'
            if not os.path.exists(outfile):
                f.write(file + '\n')
