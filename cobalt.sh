set -eo pipefail

Tbam=$1
Nbam=$2
Ts=$3
Ns=$4
java -Xmx8G -Djava.io.tmpdir=. -cp /share/work1/wangrr/local/simple/bin/cobalt-1.7.jar com.hartwig.hmftools.cobalt.CountBamLinesApplication \
    -reference $Ns \
    -reference_bam $Nbam \
    -tumor $Ts \
    -tumor_bam $Tbam \
    -output_dir cobalt_$Ts \
    -threads 16 \
    -gc_profile /share/Oncology/RD_HCC/IBFC2018028/WgsReSeq/Purple/GC_profile.hg19.1000bp.cnp
