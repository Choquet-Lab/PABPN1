import pandas as pd
import sys

bed_genes = pd.read_table(str(sys.argv[1]), header = None, dtype={0:'str', 8:'str'})
bed_genes.columns = ['chrom', 'start', 'end', 'read_id', 'score', 'strand', 'cigar', 'polya_length', 'chrom_gene', 'start_gene', 'end_gene', 'Geneid', 'score_gene', 'strand_gene', 'length']

bed_genes_count = bed_genes.groupby(by = ['chrom_gene', 'start_gene', 'end_gene', 'strand_gene', 'Geneid'])['read_id'].count().reset_index()

bed_genes_count.to_csv(str(sys.argv[2]), sep = '\t')
