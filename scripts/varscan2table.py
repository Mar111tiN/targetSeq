from os import system as shell
import os
w = snakemake.wildcards
config = snakemake.config


input = snakemake.input.vcf
output = snakemake.output
params = snakemake.params
vcf2table = params.vcf2table
varscan2table = params.varscan2table
log = snakemake.log
build = config['ref']['build']
ref = os.path.join(config['paths']['mystatic'], config['ref'][build]['genome'])


# standard output
print('input:', input)
out = input.replace('vcf', 'table')
out_ln = input.replace('vcf', 'ln.vcf')
# for indel realignment, files have to be bgzipped and indexed with tabix
shell(f"bgzip < {input} > {input}.gz")
shell(f"tabix {input}.gz")
shell(f"bcftools norm -f {ref} -o {out_ln} {input}.gz")
shell(f"cat {out_ln} | {vcf2table} > {out}")
# cleanup
shell(f"rm {input}.gz; rm {input}.gz.tbi; rm {input}; mv {out_ln} {input}")


# CONCAT THE FILES

########### !! does not work for empty files !! ### --> workaround:
# no header and merge files with fixed varscan header in merge_table


# # get first line
# head_cmd = f"mawk 'NR == 1 {{print}}' < {output_files[0]} > {output}"
# print(head_cmd)
# shell(head_cmd)
# # concat and sort the files and append to output
# cat_cmd = f"cat {' '.join(output_files)} | sort -V -k1,2 | mawk 'NR > 2 {{ print }}' >>  {output}"
# print(cat_cmd)
# shell(cat_cmd)
####################################################################

# concat and sort the files and write to output
cat_cmd = f"cat {out} | sort -V -k1,2  > {output}" #| mawk 'NR > 2 {{ print }}'
print(cat_cmd)
shell(cat_cmd)

shell(f"rm -f {out}")
print(f"Edited {input} to {output} for annotation.")
# shell(f"rm {' '.join(output_files)}")
