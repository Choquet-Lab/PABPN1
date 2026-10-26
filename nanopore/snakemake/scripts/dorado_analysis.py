import pandas as pd
import numpy as np
import re
import sys

data = open(str(sys.argv[1]), 'r')
lines = data.readlines()
reads = []
polya_list = []
for line in lines :
    try :
        polya_list.append(line.split('pt:i:')[1].split('\t')[0].strip('\n'))
        reads.append(line.split('\t')[0])
    except :
        print(line)


df = pd.DataFrame()
df['reads'] = reads
df['polya'] = polya_list
df.to_csv(str(sys.argv[2]), sep='\t')
