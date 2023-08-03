{sequenza-utils} bam2seqz -n normal.bam -t tumor.bam -gc {pipeDir}/hg19.gc50Base.txt.gz -F genome.fa -S $(which samtools) | {sequenza-utils} seqz_binning -w 50 -s - | gzip > sample.samll_seqs.gz
Rscript {pipeDir}/sequenza.R sample.samll_seqs.gz {sample}
