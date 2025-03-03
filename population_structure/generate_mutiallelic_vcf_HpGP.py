# genrating VCF files
# msa2vcf.py
import collections

dir = "/media/nosada/ToshibaMN12/pylori/multiple_fasta/HpGPcore/"


fastalist = dir+"fasta_nuc_list"
# fastalist = "temp_list"

def define_alt (allele_list):
#     ref = allele_list[0]
    alt_list = []
#     alt_list.append(ref)
    for i in range (1,len(allele_list)):
        if allele_list[i] != 'N' and allele_list[i] != '-' and allele_list[i] != ref and allele_list[i] not in alt_list and allele_list[i] != 'X':
            alt_list.append(allele_list[i])
    return alt_list
    
### generationg an alignement index
#### align nucleotide fasta using amino acid fasta
from Bio import SeqIO
from Bio.Seq import Seq
import re

strain_list = []
for strain in open("io/strain_list_HpGP_nondup_sorted.txt"):
    strain_list.append(strain[:-1])
    
vcfhead = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
header = '\t'.join(vcfhead)
samples = '\t'.join(strain_list)

print(len(strain_list))



out_vcf = open("../script/io/pylori_HpGP_decomposed.vcf", "w")
# out_vcf_allsites = open(dir+"../plink/pylori_D124_tblastn_allsites.vcf", "w")

# out_annotation = open(dir+"../io/pylori_D124_tblastn.site.txt", "w")
# out_concat_fasta = open(dir+"../plink/merged.fasta", "w")

out_vcf.write(header+"\t"+samples+"\n")

site_list = []

missing = open(dir+"missing.txt", "w")

concat_fasta = {}
position_id = 1
pos_sum = 0
### temp.txt is a list of fasta files per gene
for line in open(fastalist):
    file = line[:-1]
    print(file)
    basename = file.replace('.nuc.aligned.fasta', '')
    print(basename)
    HP = basename
    HP = HP.replace('../multiple_fasta/HpGPcore/', '')
    basename = basename.replace('../multiple_fasta/HpGPcore/', '')
    nuc_file =  dir+basename+".nuc.aligned.fasta"
    nuc_records = list(SeqIO.parse(nuc_file, 'fasta'))
    
    nuc_seq_dict = {}

    for index, nuc_record in enumerate(nuc_records):
        nuc_seq_dict[nuc_record.id] = str(nuc_record.seq.upper())

        
    offset = 0
    ref_sample = ''
    if 'HP26695' in nuc_seq_dict:
        ref_seq = nuc_seq_dict['HP26695']
        ref_sample = 'HP26695'
    else:
        missing.write(basename+" \n")
        continue
    for i in range(len(ref_seq)):
        allele_list = []
        pos_sum += i
        nuc_list = []
        
#         nuc_list.append(nuc_seq_dict[ref_sample][i])
        ref_allele = nuc_seq_dict[ref_sample][i]
   
        nX = 0
        filtered = 0
        for nuc in nuc_list:
#             print(nuc)
            if nuc == 'X' or nuc == '-':
                nX += 1
#         print(pos, nX, len(nuc_list))
        if nX >= 0.5 * len(nuc_list):# filter if more than 50% of nuculeotides are X or gaps
            filtered = 1
#         print(len(strain_list), len(nuc_list))
#         print(ref_allele, alt, nuc_list)
        pos = i - offset +1
        gappos = i + 1
         
        if ref_allele == '-':
            ID = ref_sample+':'+HP+':-:-'
            # out_annotation.write(ID+"\t"+HP+":"+str(gappos)+"\t"+str(filtered)+"\n")
            offset += 1
            continue
        for strain in strain_list:
            if strain in nuc_seq_dict:
                nuc_list.append(nuc_seq_dict[strain][i])
            else:
                nuc_list.append('N')
#         allele_list = define_alt(nuc_list)
        nuc_list_clean = []
        for nuc in nuc_list:
            if nuc == 'A' or nuc == 'T' or nuc == 'C' or nuc == 'G':
                nuc_list_clean.append(nuc)
        c = collections.Counter(nuc_list_clean)
        if len(c) <= 1:
            continue
#         print(HP, pos, c)
        allele_list, counts = zip(*c.most_common())
        ref_allele = allele_list[0]
        alt = allele_list[1:]
#         nuc_list = nuc_list[1:]
        ref_string = allele_list[0]
        if alt:### skip nonvariant sites
            for i in range(len(alt)):
                nuc_list_bin = nuc_list.copy()
                alt_string = alt[i]
                for j in range(len(nuc_list)):
                    if nuc_list_bin[j] == alt_string:
                        nuc_list_bin[j] = '1'
                    elif nuc_list[j] == ref_allele:
                        nuc_list_bin[j] = '0'
                    else:
                        nuc_list_bin[j] = '.'

                GT = '\t'.join(nuc_list_bin)
                ID = ref_sample+':'+HP+':'+str(pos)+":"+ref_allele+alt_string
                out_vcf.write("chr1\t"+str(position_id)+"\t"+ID+"\t"+ref_allele+"\t"+alt_string+'\tQUAL\tPASS\tINFO\tGT\t'+\
                          GT+"\n")
                position_id += 1

out_vcf.close() 
missing.close()
# out_concat_fasta.close()
# out_annotation.close()

