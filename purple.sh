set -eo pipefail
Tbam=$1
Nbam=$2
Ts=$3
Ns=$4

source /share/work1/wangrr/local/miniconda3/envs/initialize 

shelldir=/share/Oncology/RD_HCC/IBFC2018028/WgsReSeq/Purple/
sh $shelldir/amber.sh $Tbam $Nbam $Ts $Ns &
sh $shelldir/cobalt.sh $Tbam $Nbam $Ts $Ns &
wait
#
java -jar -Djava.io.tmpdir=. /share/work1/wangrr/local/simple/bin/purple-2.34.jar \
   -reference $Ns \
   -tumor $Ts \
   -output_dir purple_$Ts \
   -amber amber_$Ts \
   -cobalt cobalt_$Ts \
   -gc_profile /share/Oncology/RD_HCC/IBFC2018028/WgsReSeq/Purple/GC_profile.hg19.1000bp.cnp \
   -ref_genome /share/work1/wangrr/DB/hg19/hg19AddVirus.fa \
   -somatic_vcf /share/Oncology/RD_HCC/IBFC2018028/WgsReSeq/Purple/vcf/$Ts.overlap.vcf \
   -circos /share/public/software/circos/bin/circos
