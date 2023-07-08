set -eo pipefail

Tbam=$1
Nbam=$2
Ts=$3
Ns=$4
java -Xmx32G -Djava.io.tmpdir=. -cp amber-3.0.jar com.hartwig.hmftools.amber.AmberApplication \
    -reference $Ns \
    -reference_bam $Nbam \
    -tumor $Ts \
    -tumor_bam $Tbam \
    -output_dir amber_$Ts \
   -threads 16 \
   -loci GermlineHetPon.hg19.vcf
