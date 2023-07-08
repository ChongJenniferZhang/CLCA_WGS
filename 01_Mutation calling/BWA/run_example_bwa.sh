{bwa} mem -t 4 -R "@RG\tID:{sample}_{lane}\tPL:illumina\tSM:{sample}" {genome} {r1} {r2} | {samtools} view -bS - | {samtools} sort - -m 4G -T {prx}.tmp -o {prx}.sort.bam
