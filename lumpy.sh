sample=$1
Tbam=$2
Nbam=$3

smoove call -x --name $sample --exclude btu356_LCR-hs37d5.bed --fasta hg19.fa -p 8 --genotype $Tbam $Nbam

perl lumpySomatic.pl -i $sample-smoove.genotyped.vcf.gz > $sample.vcf
