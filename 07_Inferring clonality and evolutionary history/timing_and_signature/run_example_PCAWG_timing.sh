set -eo pipefail

for vcf in *.vcf
do
	sample=$(basename $vcf | sed 's/.vcf//')
	if [ ! -s out/$sample\_timed_segments.txt ]
	then
Rscript PCAWG_timing.R $vcf ./ ./
fi
done
