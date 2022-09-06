import re
import typer
import numpy as np
import pandas as pd
import pyranges as pr


def load_faidx(faidx):
    chrom2len = {}
    with open(faidx) as f:
        for line in f:
            tline = line.strip().split()
            chrom2len[tline[0]] = int(tline[1])
    return chrom2len

def load_gff(gfffile):
    re_geneid = re.compile('gene_id=([^;]+)')
    df_gff = []
    with open(gfffile) as f:
        for line in f:
            if line[0] != '#':
                tline = line.strip().split()
                chrom = tline[0]
                start = int(tline[3]) - 1 # convert to 0-base
                end = int(tline[4]) # convert to 0-base
                strand = tline[6]
                feature = tline[2]
                geneid = re_geneid.search(tline[8]).group(1)
                df_gff.append([chrom, start, end, strand, feature, geneid])
    df_gff = pd.DataFrame(df_gff,
                          columns=['Chromosome', 'Start', 'End', 'strand', 'feature', 'geneid'])
    return df_gff

def main(gfffile: str = typer.Argument(..., help="gff3注释文件"),
         faidx: str = typer.Argument(..., help="samtools faidx文件"),
         updownstream_length: int = typer.Argument(2000, help="上下游的范围"),
         outfile: str = typer.Argument(..., help="输出文件名")):
    chrom2len = load_faidx(faidx)
    df_gff = load_gff(gfffile)
    
    # 汇总每个基因上各个feature的长度
    features = {'exon', 'intron', 'CDS', 'three_prime_UTR', 'five_prime_UTR'}
    fdf = (df_gff.query('feature in @features')[['Chromosome', 'Start', 'End', 'geneid', 'feature']]
                 .assign(length=lambda x: x['End'] - x['Start'])
                 .drop(columns=['Chromosome', 'Start', 'End'])
                 .groupby(['geneid', 'feature']).sum().reset_index()
                 .pivot(index='geneid', columns='feature', values='length')
                 .reset_index()
                 .fillna(0))
    # 根据annovar注释优先级：
    fdf['exon'] = fdf['exon'] - fdf[['three_prime_UTR', 'five_prime_UTR']].sum(axis=1)
    
    # 计算上下游范围
    # 正链
    gdf_pos = df_gff.query('(feature == "gene") and (strand == "+")').copy()
    gdf_pos['upstream_start'] = gdf_pos['Start'] - updownstream_length
    gdf_pos['upstream_end'] = gdf_pos['Start']
    gdf_pos['downstream_start'] = gdf_pos['End']
    gdf_pos['downstream_end'] = gdf_pos['End'] + updownstream_length
    # 负链
    gdf_neg = df_gff.query('(feature == "gene") and (strand == "-")').copy()
    gdf_neg['upstream_start'] = gdf_neg['End']
    gdf_neg['upstream_end'] = gdf_neg['End'] + updownstream_length
    gdf_neg['downstream_start'] = gdf_neg['Start'] - updownstream_length
    gdf_neg['downstream_end'] = gdf_neg['Start']
    
    # 合并正负链区间
    cols_upstream = ['Chromosome', 'upstream_start', 'upstream_end', 'geneid']
    rename_upstream = {'upstream_start': 'Start', 'upstream_end': 'End'}
    up_pos = gdf_pos[cols_upstream].rename(columns=rename_upstream).assign(feature='upstream')
    up_neg = (gdf_neg[cols_upstream].rename(columns=rename_upstream).assign(feature='upstream'))

    cols_downsteam = ['Chromosome', 'downstream_start', 'downstream_end', 'geneid']
    rename_downstream = {'downstream_start': 'Start', 'downstream_end': 'End'}
    down_pos = (gdf_pos[cols_downsteam].rename(columns=rename_downstream).assign(feature='downstream'))
    down_neg = (gdf_neg[cols_downsteam].rename(columns=rename_downstream).assign(feature='downstream'))

    updown_gdf = pd.concat([up_pos, up_neg, down_pos, down_neg], ignore_index=True).sort_values(['Chromosome', 'Start'])
    print()
    # 处理超过染色体边界的
    updown_gdf['Start'] = updown_gdf['Start'].map(lambda x: x if x >= 0 else 0)
    
    chrom_edges = updown_gdf['Chromosome'].map(chrom2len).values + 1
    selection_outrange = updown_gdf['End'].values > chrom_edges
    updown_gdf.loc[selection_outrange, 'End'] = chrom_edges[selection_outrange]
    
    updown_gdf['Length'] = updown_gdf['End'].values - updown_gdf['Start'].values
    length_upstream = updown_gdf.query('feature == "upstream"').set_index('geneid')['Length']
    length_downstream = updown_gdf.query('feature == "downstream"').set_index('geneid')['Length']
    
    # 计算与基因相交的区间
    gr_updownstream = pr.pyranges.PyRanges(updown_gdf)
    gr_genes = pr.pyranges.PyRanges(df_gff.query('feature == "gene"')[['Chromosome', 'Start', 'End']])
    overlap_upstream = {}
    overlap_downstream = {}
    for (geneid, feature), idf in (gr_updownstream.intersect(gr_genes, strandedness=False)
                                      .cluster(by=['geneid', 'feature'], strand=False)
                                      .df.groupby(['geneid', 'feature'], sort=False)):
        if feature == 'upstream':
            overlap_upstream[geneid] = idf['End'].max() - idf['Start'].min()
        elif feature == 'downstream':
            overlap_downstream[geneid] = idf['End'].max() - idf['Start'].min()
    
    # 添加上下游长度
    fdf['upstream'] = fdf['geneid'].map(lambda x: length_upstream.get(x, 0) - overlap_upstream.get(x, 0))
    fdf['downstream'] = fdf['geneid'].map(lambda x: length_downstream.get(x, 0) - overlap_downstream.get(x, 0))
    fdf.to_csv(outfile, index=False, sep='\t')
    

if __name__ == '__main__':
    typer.run(main)