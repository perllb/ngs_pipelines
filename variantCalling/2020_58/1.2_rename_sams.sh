

cd /data/bnf/dev/per/proj/TWIST/2020_58_Hedenfalk/bam

#
# Rename files
#
for file in *.sam; do echo $file; newname=$(echo CTG-${file}); echo $newname; mv $file $newname; done

