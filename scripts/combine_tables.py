import pandas as pd


# ############## SNAKEMAKE ################################
w = snakemake.wildcards
config = snakemake.config

input = snakemake.input
output = snakemake.output
threads = snakemake.threads
log = snakemake.log


# ############## ANNOVAR FILE ########################################
#####################################################################


print(f'Loading annovar file and adjusting columns {input.anno_edit}.')
anno_df = pd.read_csv(input.anno_edit, sep='\t', low_memory=False, dtype={'Chr': str, 'Start': int, 'End': int}, na_values=['---', 'NaN', 'NAN']).fillna('.').sort_values(['Chr', 'Start'])

anno_df['TVAF'] = (anno_df['TR2'] / (anno_df['TR1'] + anno_df['TR2'])).round(3)
anno_df = anno_df.rename(columns={'readDepth': 'Tdepth'})
# RENAMING ##################################

# get rename dict for .refGene annotation
refgen_dict = {col: col.replace(".refGene", "") for col in anno_df.columns if ".refGene" in col}
# rename the columns
anno_df = anno_df.rename(columns=refgen_dict)


# resort the columns
cols = list(anno_df.columns)
# for i, col in enumerate(cols):
#   print(i, '-', col)

start_cols = cols[:16]
quant_cols = ['Tdepth'] + ['TVAF'] + cols[19:25]
anno_cols = cols[35:]


# MERGE Func into ExonicFunc ##################################
# merge ExonFunc from Func, if ExonFunc not available
anno_df.loc[anno_df['ExonicFunc'] == ".", 'ExonicFunc'] = anno_df['Func']


# ############ -->COLS #####################################
new_cols = start_cols + quant_cols
# ############## MERGE WITH EB AND FISHER ############################
#####################################################################

print(f'Loading FisherScore from file {input.fisher} and merging into annotation.')
fisher_df = pd.read_csv(input.fisher, sep='\t', dtype={'Chr': str, 'Start': int, 'End': int}).fillna('.').sort_values(['Chr', 'Start'])
merge_F = anno_df.merge(fisher_df, on=['Chr', 'Start', 'End', 'Ref', 'Alt'])
fisher_cols = list(fisher_df.columns[5:])
new_cols += fisher_cols     # -->COLS

print(f'Loading HDR-data from file {input.HDR} and merging into annotation.')
HDR_df = pd.read_csv(input.HDR, sep='\t')
merge_FH = pd.merge(merge_F, HDR_df, how='inner', on=['Chr', 'Start', 'End', 'Ref', 'Alt', 'Gene'])
HDR_cols = list(HDR_df.columns[6:])
new_cols += HDR_cols

print(f'Loading Primer3-data from file {input.primer} and merging into annotation.')
primer_df = pd.read_csv(input.primer, sep='\t')
merge_FHP = pd.merge(merge_FH, primer_df, how='inner', on=['Chr', 'Start', 'End', 'Ref', 'Alt', 'Gene'])
primer_cols = list(primer_df.columns[6:])
new_cols += anno_cols + primer_cols     # -->COLS


# ############## REARRANGE COLUMNS ###################################
#####################################################################

anno_combined = merge_FHP[new_cols]
anno_combined.to_csv(output[0], sep='\t', index=False)
print(f'Writing combined annotation to {output[0]}.')
# resort the columns
# cols = list(anno_combined.columns)
# for i, col in enumerate(cols):
#   print(i, '-', col)