r1=$1
r2=$2
prx=$3
mkdir -p $prx
STAR --runThreadN 4 --quantMode TranscriptomeSAM GeneCounts  --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20  --alignMatesGapMax 1000000   --twopassMode Basic --genomeDir /hg19/STARgenomeDir --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --sjdbGTFfile /RNA/hg19_E.gtf --outReadsUnmapped None --chimSegmentMin 12  --chimJunctionOverhangMin 12    --alignIntronMax 100000 --chimSegmentReadGapMax parameter 3  --alignSJstitchMismatchNmax 5 -1 5 5  --readFilesIn  $r1 $r2 --outFileNamePrefix  $prx/$prx --chimOutType Junctions

