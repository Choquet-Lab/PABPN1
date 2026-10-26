import pandas as pd
import sys
import os

input = sys.argv[1]
output = sys.argv[2]

splicing = pd.read_table(str(input))
splicing_counts = splicing.groupby(['chrom', 'intron_start', 'intron_end', 'strand', 'gene_name', 'intron_count', 'splice_status'])['read_name'].count().reset_index()
splicing_counts_piv = splicing_counts.pivot_table(index=['chrom', 'intron_start', 'intron_end', 'strand', 'gene_name', 'intron_count'], columns='splice_status',values='read_name').reset_index().fillna(0)
splicing_counts_filt = splicing_counts_piv[splicing_counts_piv['NO'] + splicing_counts_piv['YES'] > 5].reset_index(drop=True).drop(columns=['SKP', 'SKP3SS', 'SKP5SS', 'UND'])
splicing_counts_filt['fraction_spliced'] = splicing_counts_filt['YES'] / (splicing_counts_filt['YES'] + splicing_counts_filt['NO'])
splicing_counts_filt.to_csv(str(output), sep='\t')
