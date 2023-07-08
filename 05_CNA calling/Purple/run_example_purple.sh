set -eo pipefail
Tbam=T.bam
Nbam=N.bam
Ts=T
Ns=N


shelldir=$(pwd)
sh $shelldir/amber.sh $Tbam $Nbam $Ts $Ns &
sh $shelldir/cobalt.sh $Tbam $Nbam $Ts $Ns &
wait
#
java -jar -Djava.io.tmpdir=. purple-2.34.jar \
   -reference $Ns \
   -tumor $Ts \
   -output_dir purple_$Ts \
   -amber amber_$Ts \
   -cobalt cobalt_$Ts \
   -gc_profile GC_profile.hg19.1000bp.cnp \
   -ref_genome hg19AddVirus.fa 
