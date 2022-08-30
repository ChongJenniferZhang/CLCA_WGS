{sequenza-utils} bam2seqz -n normal.bam -t tumor.bam -gc {pipeDir}/hg19.gc50Base.txt.gz -F /share/public/database/hg19/hg19.fa -S /share/public/software/samtools/0.1.18/samtools | {sequenza-utils} seqz_binning -w 50 -s - | gzip > sample.samll_seqs.gz
Rscript {pipeDir}/sequenza.R sample.samll_seqs.gz {sample}
