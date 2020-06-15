import os
import re
import yaml
import argparse
import math

# ############ SETUP ##############################
configfile: "configs/config_AMLValidation.yaml"
# configfile: "configs/config.json"
workdir: config['workdir']

# include helper functions
include: "includes/io.snk"
include: "includes/utils.snk"

# FILE SEARCH ###########
#   included_files = tumor_normal_pairs = sample_types = None
#  if not included_files:  # do this only once
#  included_files, tumor_normal_pairs, sample_types = get_files(config['inputdir'])
# included files : list of (fastq_path, sample, type, read, tumor/normal) tuples
# tumor_normal_pairs : list of 'sample_tumor_normal' strings
# sample_types : list of 'sample_type' strings


# retrieve the file_df with all the file paths from the samplesheet
sheet_file = config['samplesheet']
# remove the header and other unneccessary info
shell(f"cat {sheet_file} | sed -n '/Data/,$ p' > sample_short.csv")
# make global samples variable
samples = list(pd.read_csv('sample_short.csv', sep=',', skiprows=1)['Sample_Name'])
print('Running samples:', ', '.join(samples))

chrom_list = get_chrom_list(config)


# ############ INCLUDES ##############################
include: "includes/demulti.snk"
# include: "includes/fastq.snk"
include: "includes/QC.snk"
include: "includes/map.snk"
include: "includes/processBAM.snk"
include: "includes/dedup.snk"
include: "includes/umi_filter.snk"
include: "includes/freebayes.snk"
include: "includes/annotate.snk"
include: "includes/filter.snk"
include: "includes/filterbam.snk"

# convenience variables
ref_gen = full_path('genome')
# specified wildcards have to match the regex
wildcard_constraints:
    # eg sample cannot contain _ or / to prevent ambiguous wildcards
    sample = "[^/.]+",
    # read = "[^_/.]+",
    # read_or_index = "[^_/.]+",
    # filter = "filter[0-9]+"
    # folder = "^((?!filter).)*$"


# ############## MASTER RULE ##############################################

rule all:
    input:
        expand("umi/{sample}.UF.bam", sample = samples),
        expand("filter/{sample}.filter1.csv", sample = samples),
        expand("filterbam/{sample}.filter1.txt", sample = samples),
        # "QC/libraryQC.html",
        # "QC/insertQC.html",
        "barcodes/unmatched_barcodes.txt",
        "QC/fastQC.html",
        # expand("coverBED/{sample}.txt", sample = samples)

###########################################################################


# print out of installed tools
onstart:
    print("    TWIST TARGETED SEQUENCING PIPELINE STARTING.......")

    ##########################
    # shell("echo Conda-environment: $CONDA_PREFIX")
    # shell('echo $PATH')
    # write config to the results directory
    path_to_config = os.path.join(config['workdir'], "config.yaml")
    with open(path_to_config, 'w+') as stream:
        yaml.dump(config, stream, default_flow_style=False)
    # create logs folder
#     shell("conda list | show_awk")
#     shell("ls -l ${{TOOLS}}/bulltools/links | sed -nr '/->/s_.* ([^/]+ -> .*)_  \\1_p' ")
    # create scratch/tmp/hWES folder for storing stuff


onsuccess:
    # shell("export PATH=$ORG_PATH; unset ORG_PATH")
    print("Workflow finished - everything ran smoothly")

    # cleanup
    if config['setup']['cleanup']['run']:
        no_bams = config['setup']['cleanup']['keep_only_final_bam']

        split_fastq_pattern = '\.[0-9]+\.fastq.gz'
        split_bam_pattern = 'chr[^.]+\..*'

        # remove split fastqs
        shell("rm -rf ubam realigned bam_metrics insert_metrics pileup varscan fastqc mapped bam_done")


        shell("ls fastq | grep -E '{split_fastq_pattern}' | sed 's_^_fastq/_' | xargs -r rm -f")
        
        # remove split recalib bams
        shell("ls recalib | grep -E '{split_bam_pattern}' | sed 's-^-recalib/-' | xargs -r rm -f")

        if no_bams:
            shell("rm -r bam_merge ")
        else:
            # only remove split_bams in bam_merge
            shell("ls bam_merge | grep -E '{split_bam_pattern}' | sed 's-^-bam_merge/-' | xargs -r rm -f")
