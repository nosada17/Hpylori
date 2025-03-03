# Description of how the data was processed

Alignment files were processed using `generate_multiallelic_vcf_HpGP.py` script. Output file name is `pylori_HpGP_decomposed.vcf`

The vcf file was processed using the following commands.
```
plink --vcf io/pylori_D7_decomposed.vcf \
    --mind 0.1 --geno 0.2 --maf 0.01 \
    --remove ../plink/pylori_D7_prunned_pihat_out.txt \
    --make-bed --keep-allele-order \
    --out ../plink/pylori_D7_decomposed_pihat_missing --const-fid
```
