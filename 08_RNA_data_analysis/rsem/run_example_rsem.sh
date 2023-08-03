bam=$1
sample=$2
mkdir -p $sample
export PATH=/RSEM-1.3.3/bin/:$PATH
cd $sample
rsem-calculate-expression --bam --no-bam-output -p 4 --paired-end  --forward-prob 0 $bam //hg19/STARgenomeDir/hg19AddVirus $sample

