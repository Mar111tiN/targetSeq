import os
import pandas as pd
from scipy.stats import fisher_exact as fe
import math
import numpy as np
from multiprocessing import Pool
from script_utils import show_output, show_command


w = snakemake.wildcards
config = snakemake.config
input = snakemake.input
output = snakemake.output
threads = snakemake.threads
log = snakemake.log


#####################################################################
# ###################### FISHER SCORE ################################

def get_fisher_exact(row):
    T1plus = row['TR1+']
    T1min = row['TR1'] - T1plus
    T2plus = row['TR2+']
    T2min = row['TR2'] - T2plus
    mat = np.matrix([[T1plus, T2plus], [T1min, T2min]])
    fisher_p = fe(mat)[1]
    if fisher_p:
        return round(-10*math.log(fisher_p, 10), 1)
    return 5000

print(f'Computing FisherScore for file {input[0]}')
df = pd.read_csv(input[0], sep='\t')
cols = list(df.columns)[:5] 
df['FisherScore'] = df.apply(get_fisher_exact, axis=1)

# reduce to important cols
cols += ['FisherScore']
# write file to filtered
df[cols].to_csv(output[0], sep='\t', index=False)
show_output(f'FisherScore for file {input[0]} written to {output[0]}', color='success')