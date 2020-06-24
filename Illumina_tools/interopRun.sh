#!/bin/sh

source /sw/easybuild/software/Core/Miniconda3/4.5.12/etc/profile.d/conda.sh

conda activate /home/per/.conda/envs/Illumina_interop

python /fs1/per/github/ngs_pipelines/Illumina_tools/Interop_STATS.py $1 > $1/Interop_STATS.txt

cat $1/Interop_Stats_*.csv > $1/Interop_STATS.All.csv

