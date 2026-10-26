import pandas as pd
import scipy
import os
import sys
from scipy import stats
from statsmodels.stats.multitest import multipletests

##############################################
##############################################
##############################################
##############################################

compartment = sys.argv[1]
list_splicing = os.listdir(sys.argv[2])
gene_name_ref = sys.argv[3]
intron_ref = sys.argv[4]

list_shctl = []
list_shpab = []

for file in list_splicing :
    if compartment in file :
        if ('shCTL' in file) or ('shNT' in file) :
            list_shctl.append(pd.read_table(sys.argv[2] + '/' + file, dtype={'chrom':'str'}))
        elif ('shPABPN1' in file) :
            list_shpab.append(pd.read_table(sys.argv[2] + '/' + file, dtype={'chrom':'str'}))

shCTL = pd.concat(list_shctl).reset_index(drop=True)
shPAB = pd.concat(list_shpab).reset_index(drop=True)


def compare_shCTL_to_shPAB_splicing(shCTL_df, shPAB_df, comp):
    # Count number of reads per splicing status per intron (cyto)
    shCTL_counts = pd.DataFrame(shCTL_df.groupby(['chrom', 'intron_start', 'intron_end', 'strand', 'gene_name', 'intron_count', 'splice_status'])['read_name'].count()).reset_index()
    shCTL_counts_piv = shCTL_counts.pivot_table(index=['chrom', 'intron_start', 'intron_end', 'strand', 'gene_name', 'intron_count'], columns='splice_status',values='read_name').reset_index().fillna(0)

    shPAB_counts = pd.DataFrame(shPAB_df.groupby(['chrom', 'intron_start', 'intron_end', 'strand', 'gene_name', 'intron_count', 'splice_status'])['read_name'].count()).reset_index()
    shPAB_counts_piv = shPAB_counts.pivot_table(index=['chrom', 'intron_start', 'intron_end', 'strand', 'gene_name', 'intron_count'], columns='splice_status',values='read_name').reset_index().fillna(0)

    # Consider only YES and NO for now
    shCTL_counts_filt = shCTL_counts_piv[shCTL_counts_piv['NO'] + shCTL_counts_piv['YES'] > 20].reset_index(drop=True).drop(columns=['SKP', 'SKP3SS', 'SKP5SS', 'UND'])
    shPAB_counts_filt = shPAB_counts_piv[shPAB_counts_piv['NO'] + shPAB_counts_piv['YES'] > 20].reset_index(drop=True).drop(columns=['SKP', 'SKP3SS', 'SKP5SS', 'UND'])

    # Merge the counts from the two conditions
    counts = pd.merge(shCTL_counts_filt, shPAB_counts_filt,on=['chrom', 'intron_start', 'intron_end', 'strand', 'gene_name', 'intron_count'])

    if comp == "cyto":
        # Keep only introns that have some unspliced reads in at least one condition
        counts_filt = counts[(counts['NO_x'] > 3) | (counts['NO_y'] > 3)].reset_index(drop=True).rename(columns={'gene_name': 'gene'})
    elif comp == "noyau":
        # Keep only introns that have both spliced and unspliced reads in at least one condition
        counts_filt = counts[((counts['YES_x'] >= 3) | (counts['YES_y'] >= 3)) & ((counts['NO_x'] >= 3) | (counts['NO_y'] >= 3))].reset_index(drop=True).rename(columns={'gene_name': 'gene'})

    # Add gene names
    counts_filt['gene_id'] = counts_filt['gene'].str.split("\\.").str[0]
    gene_names_df = pd.read_table(gene_name_ref, header = None)
    gene_names_df.columns = ['gene_name', 'gene_id']
    counts_filt = counts_filt.merge(gene_names_df, on='gene_id')

    # Get intron annotations
    intron_df = pd.read_table(intron_ref, dtype={"chrom": "string"})
    #intron_df['chrom'] = intron_df['chrom'].str[3:]

    intron_fields = ['chrom', 'start', 'end', 'strand', 'count', 'gene', 'intron_order']
    counts_filt = counts_filt.merge(intron_df,
                                    left_on=['chrom', 'intron_start', 'intron_end', 'strand', 'gene', 'intron_count'],
                                    right_on=['chrom', 'start', 'end', 'strand', 'gene', 'count'])

    # Perform statistical testing
    counts_filt['pvalue'] = counts_filt.apply(lambda row: scipy.stats.fisher_exact([[row.NO_x, row.YES_x], [row.NO_y, row.YES_y]])[1], axis=1)
    counts_filt = counts_filt.sort_values(by='pvalue').drop(columns=['gene', 'gene_id', 'count']).drop_duplicates().reset_index(drop=True).fillna('single')
    counts_filt['intron_order'] = pd.Categorical(counts_filt['intron_order'],categories=['last', 'single', 'middle', 'first'], ordered=True)
    counts_filt_bis = counts_filt.sort_values('intron_order').drop_duplicates(subset=['chrom', 'intron_start', 'intron_end', 'strand']).sort_values(by='pvalue').reset_index(drop=True)
    counts_filt_bis['FDR'] = multipletests(counts_filt_bis['pvalue'].astype(float), alpha=0.05, method='fdr_bh')[1]

    # Measure fraction spliced
    counts_filt_bis['fraction_spliced_x'] = counts_filt_bis['YES_x'] / (counts_filt_bis['YES_x'] + counts_filt_bis['NO_x'])
    counts_filt_bis['fraction_spliced_y'] = counts_filt_bis['YES_y'] / (counts_filt_bis['YES_y'] + counts_filt_bis['NO_y'])
    counts_filt_bis['delta'] = counts_filt_bis['fraction_spliced_x'] - counts_filt_bis['fraction_spliced_y']

    counts_filt_bis['dir'] = 'no_diff'
    counts_filt_bis.loc[counts_filt_bis['delta'] > 0, 'dir'] = 'more_IR_PABPN1'
    counts_filt_bis.loc[counts_filt_bis['delta'] < 0, 'dir'] = 'less_IR_PABPN1'

    counts_filt_bis['sig'] = 'no'
    counts_filt_bis.loc[(counts_filt_bis['FDR'] < 0.05) & (abs(counts_filt_bis['delta']) > 0.05), 'sig'] = 'yes'

    return (counts_filt_bis)

counts_filt_bis = compare_shCTL_to_shPAB_splicing(shCTL,shPAB,compartment)

counts_filt_bis.to_csv(sys.argv[5], sep = '\t')
