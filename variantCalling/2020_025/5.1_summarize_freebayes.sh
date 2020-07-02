#!/bin/sh

baseOut="/data/bnf/dev/per/proj/TWIST/2020_025_JKarlsson"
scrdir="/trannel/dev/per/proj/TWIST/2020_025_JKarlsson/scr"
FBdir=${baseOut}/vcf_freebayes_cohort

# create xls

## freebayes
/data/bnf/sw/vcflib/bin/vcfbreakmulti ${FBdir}/freebayes_cohortall.vcf > ${FBdir}/freebayes_cohortall_split.vcf
${scrdir}/freebayes_to_table.pl ${FBdir}/freebayes_cohortall_split.vcf > ${FBdir}/freebayes_cohortall_split.xls

