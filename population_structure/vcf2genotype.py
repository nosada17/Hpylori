# this scirpt convert vcf file to genotype file used for PCAUMAP_blobal.py

import pandas as pd
import re
import statistics
import re

dir = "directory_for_opening"

##### LD prunned vcf data, diploid format, genotypes are either 0/0, 1/1, or ./.
vcf_open = "HpGP_decomposed_pihat_missing_prunned.vcf"

with open(vcf_open) as vcf:
    line = vcf.readlines()[6]### depends on how the vcf file was created
    itemList = line[:-1].split()
    sample_names = itemList[9:]
vcf.close()

print(sample_names)

### PLINK add　"0_" to the strain names. removing them. 
for index, sample in enumerate(sample_names):
    sample_names[index] = re.sub('^0_', '', sample)
    
    
nsample = len(sample_names)
print(sample_names)

df =  pd.DataFrame(data=None, index=sample_names, columns=None, dtype=None, copy=False)

# print(df)qqq

missing = 0

#### output file
table = open("io/genotype_HpGP_decomposed_pihat_missing_prunned.txt", "w")

sample_name_string = '\t'.join(sample_names)
table.write("ID\t"+sample_name_string+"\n")

nline = 0
with open(vcf_open) as vcf:
    for line in vcf:
        nline += 1
        if re.match('\#', line):
            continue
        itemList = line[:-1].split()
        GT_list = itemList[9:]
        for index, GT in enumerate(GT_list):
            if GT == '0/0':
                GT_list[index] = '0'
            elif GT == '1/1':
                GT_list[index] = '1'
            elif GT == './.':
                GT_list[index] = '.'
        alt = itemList[4].split(',')
        try:
            mode = statistics.mode(GT_list)
        except:
            mode = 0
        if mode == '.':
            missing += 1
            continue
        for i, nuc in enumerate(alt):
            data = ''
            for GT in GT_list:
                if GT == '.':
                    data += '\t'+str(mode)
                elif int(GT) == i+1:
                    data += '\t1'
                else:
                    data += '\t0'
            colname = itemList[2]+"_"+nuc
            table.write(colname+data+"\n")
        if nline%10000 == 0:
            print(nline)
vcf.close()
table.close()
    
print(missing)
