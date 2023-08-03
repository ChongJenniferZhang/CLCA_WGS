gatk --java-options "-Xmx8g" Mutect2 -R hg19.fa -I T.bam -tumor Tsample -I N.bam -normal Nsample -O gatk4.vcf.gz -pon Wgs.PoN.vcf.gz -L chr --native-pair-hmm-threads 8
gatk --java-options "-Xmx2g" FilterMutectCalls -V gatk4.vcf.gz -O raw.vcf.gz -XL wgEncodeDacMapabilityConsensusExcludable.bed
