import pandas as pd
import os
from HDR_utils import masterHDR
from script_utils import show_output


# ############## SNAKEMAKE ################################
w = snakemake.wildcards
config = snakemake.config
params = snakemake.params
MINSIM = params.min_sim
min_q = params.min_q
PAD = min(config['HDR']['padding'], config['filter_bam']['padding'])
HDRMINCOUNT = params.min_HDR_count

i = snakemake.input

tumor_bam = i.bam
filter_file = i.filter_file
filter_pileup = i.pileup
out_file = snakemake.output.HDR


show_output(
    f'Starting HDR analysis of {filter_file}. [MIN_SIM={MINSIM}, PAD={PAD}]')
# GET THE mutation file for masterHDR
show_output(f'Importing {filter_file} for HDR detection', time=False)
filter_df = pd.read_csv(filter_file, sep='\t').loc[:, [
    'Chr', 'Start', 'End', 'Ref', 'Alt', 'Gene']]

HDR_df = masterHDR(
    pileup_file=filter_pileup,
    tumor_bam=tumor_bam,
    filter_df=filter_df,
    MINSIM=MINSIM,
    padding=PAD,
    min_q=min_q,
    min_HDR_count=HDRMINCOUNT
)

HDR_len = len(HDR_df.query('HDRcount > 0').index)

show_output(
    f"Found {HDR_len} possible HDR mutations. Writing to file {out_file}")
HDR_df.to_csv(out_file, sep='\t', index=False)
