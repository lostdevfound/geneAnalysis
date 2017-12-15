#!/bin/bash
###### assemble reads using SOAPdenovo ######
set -e

currentDir=`pwd`
refFasta=${currentDir}/genome-ref/chloroplast-reference-genome.fasta

workdata=( "wila" )
# workdata=( "casc" "cent" "colb" "fugg" "gold" "mtho" "nugg" "sora" "wila" )

for dirName in ${workdata[@]}
do
    seqPref=`echo $dirName | tr [:lower:] [:upper:]`
    cd workdata/${dirName}/contigs/soap

    echo "######### start assembling with SOAPdenovo $dirName #########"
    SOAPdenovo-63mer all -s ${dirName}-soap.config -K 63 -R -o ${dirName} &&

    echo "######### start aligning $dirName #########"
    bwa mem -R '@RG\tID:HWI-ST765\tSM:sample1\tLB:None\tPL:Illumina' $refFasta ${dirName}.contig > soap-${dirName}-contigs.sam

    echo "########################### sam processing start: ${dirName} ###############################"
    samtools view -F 4 -bh soap-${dirName}-contigs.sam | samtools sort -n - | samtools fixmate -m - soap-${dirName}-contigs.fix.bam &&
    samtools sort soap-${dirName}-contigs.fix.bam | samtools markdup -r - soap-${dirName}-contigs.sorted.bam &&
    samtools index soap-${dirName}-contigs.sorted.bam

    echo "######### start assembling with NOVOplasty $dirName #########"
    cd ../novoplasty

    perl ~/bin/NOVOPlasty/NOVOPlasty2.6.3.pl -c ${dirName}-novop-config.txt

    echo "######### start aligning $dirName #########"
    bwa mem -R '@RG\tID:HWI-ST765\tSM:sample1\tLB:None\tPL:Illumina' $refFasta Contigs_1_${dirName}.fasta > nvp-${dirName}-contigs.sam

    echo "########################### sam processing start: ${dirName} ###############################"
    samtools view -F 4 -bh nvp-${dirName}-contigs.sam | samtools sort -n - | samtools fixmate -m - nvp-${dirName}-contigs.fix.bam &&
    samtools sort nvp-${dirName}-contigs.fix.bam | samtools markdup -r - nvp-${dirName}-contigs.sorted.bam &&
    samtools index nvp-${dirName}-contigs.sorted.bam

    cd $currentDir
done
