#!/bin/sh

baseOut="/data/bnf/dev/per/proj/TWIST/2020_025_JKarlsson/"
vardir=${baseOut}/vcf_vardict
TNdir=${baseOut}/vcf_tnscope
FBdir=${baseOut}/vcf_freebayes_cohort

# create xls

## freebayes
/data/bnf/sw/vcflib/bin/vcfbreakmulti ${FBdir}/freebayes_cohortall.vcf > ${FBdir}/freebayes_cohortall_split.vcf
${FBdir}/freebayes_to_table.pl ${FBdir}/freebayes_cohortall_split.vcf > ${FBdir}/freebayes_cohortall_split.xls

## vardict
for vcf in ${vardir}/*vcf
do echo $vcf

   vcfbase=$(basename $vcf .vcf)

   /data/bnf/sw/vcflib/bin/vcfbreakmulti ${vardir}/${vcf} > ${vardir}/${vcfbase}_split.vcf
   ${vardir}/vardict_to_table.pl ${vardir}/${vcfbase}_split.vcf > ${vardir}/${vcfbase}_split.xls

done




