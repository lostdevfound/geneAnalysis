#!/bin/bash
if [ $1 == 'bcf' ]; then
    bcftools view $2 | zgrep -v "^##" | awk 'BEGIN{OFS="\t"} {print  $2, $4, $5, $6}' | less
elif [ $1 == 'bam' ]; then
    samtools view $2 | awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4, $5, $6}' | less
else
    echo "wrong argument try:"
    echo "usage: viewfields bcf input.bcf"
    echo "usage: viewfields bam input.bam"
fi
