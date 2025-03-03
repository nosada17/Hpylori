# Description of how the data was processed

Alignment files were processed using `generate_multiallelic_vcf_HpGP.py` script. Output file name is `pylori_HpGP_decomposed.vcf`

The vcf file was processed using the following commands.
```
# convert to plink binary format
plink --vcf pylori_HpGP_decomposed.vcf \
    --mind 0.1 --geno 0.2 --maf 0.01 \
    --make-bed --keep-allele-order \
    --out pylori_HpGP_decomposed_pihat_missing --const-fid
```
```
# LD prunning
plink --bfile pylori_HpGP_decomposed_pihat_missing --indep-pairwise 20 1 0.5
plink --bfile pylori_HpGP_decomposed_pihat_missing \
    --extract plink.prune.in --make-bed  \
    --out pylori_HpGP_decomposed_pihat_missing_prunned
```
```
# converting back to vcf format
plink --bfile ../plink/pylori_HpGP_decomposed_pihat_missing_prunned \
    --recode vcf \
    --out pylori_HpGP_decomposed_pihat_missing_prunned
```

The file `pylori_HpGP_decomposed_pihat_missing_prunned` is used for the input of `vcf2genotype.py`, which generates the input file for `PCAUMAP_global.py`
