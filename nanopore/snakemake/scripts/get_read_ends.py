import sys
import numpy as np
import pandas as pd
import pysam
from collections import Counter

import re
import math

import pybedtools
from pybedtools import BedTool

# CONFIG

pybedtools.helpers.set_bedtools_path(path='/home/sbache/projects/def-kchoquet/sbache/scripts/bedtools2/bin')

bamFile = pybedtools.BedTool(sys.argv[1])
readmap = pd.read_table(sys.argv[2], header = None)

readmap = readmap[[3,7,11,8,9,10,13]]
readmap.columns = ['reads', 'polya_length', 'Geneid', 'chrom_gene', 'start_gene', 'end_gene', 'strand_gene']

def get_read_end_bedtool(bamFile):
    bedFile = bamFile.bam_to_bed()
    bedFile_df = bedFile.to_dataframe(low_memory=False)

    read_end = []

    for i in range(0, len(bedFile_df)):

        chrom = str(bedFile_df['chrom'].iloc[i])
        start = bedFile_df['start'].iloc[i]
        end = bedFile_df['end'].iloc[i]
        read = bedFile_df['name'].iloc[i]
        score = bedFile_df['score'].iloc[i]
        strand = bedFile_df['strand'].iloc[i]

        if (strand == "-"):
            pos_1 = start
            pos_2 = start + 1

        if (strand == "+"):
            pos_1 = end - 1
            pos_2 = end

        read_end.append([chrom, pos_1, pos_2, read, score, strand])

    read_end_df = pd.DataFrame(read_end)
    read_end_df.columns = ['chrom','start','end','reads','score','strand']
    #bedFile_df.to_csv('/home/sbache/projects/def-kchoquet/sbache/test_bedfile_df.txt', sep = '\t')
    return read_end_df


# Apply function
read_ends_df = get_read_end_bedtool(bamFile)
read_ends_df = read_ends_df.merge(readmap, on='reads')

read_ends_df.to_csv(sys.argv[3], sep="\t", index=False, header=False)
