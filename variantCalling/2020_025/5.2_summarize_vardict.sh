#!/bin/sh

baseOut="/data/bnf/dev/per/proj/TWIST/2020_025_JKarlsson/"
vardir=${baseOut}/vcf_vardict

## vardict
for vcf in ${vardir}/*vcf
do echo $vcf

   vcfbase=$(basename $vcf .vcf)

   ${vardir}/vardict_to_table.pl ${vardir}/${vcfbase}_split.vcf > ${vardir}/${vcfbase}_split.xls

done





