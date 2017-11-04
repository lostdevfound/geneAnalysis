curdir=`pwd`

# dirList=( "casc" "cent" "colb" "fugg" "gold" "mtho" "nugg" "sora" "wila" )
dirList=( "fugg" )
for dir in ${dirList[@]}
do
    echo "unpacking ${dir}"
    bcftools view ${curdir}/workdata/${dir}/aln/${dir}-filterted-variance.vcf.gz > ${curdir}/workdata/${dir}/aln/${dir}-filtered-variance.vcf
done
