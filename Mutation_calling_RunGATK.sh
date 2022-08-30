{java} -Xmx5G -jar {gatk} -T RealignerTargetCreator -R {genome} -I {in} -o {prx}.intervals -known {dbsnp} -L {chr}
{java} -Xmx16G -jar {gatk} -T IndelRealigner -R {genome} -targetIntervals {prx}.intervals -I {in} -o {prx}.realn.bam -known {dbsnp} -L {chr}
{java} -Xmx16G -jar {gatk} -T BaseRecalibrator -R {genome} -I {prx}.realn.bam -knownSites {dbsnp} -o {prx}.grp -L {chr} -nct 4
{java} -Xmx16G -jar {gatk} -T PrintReads -R {genome} -I {prx}.realn.bam -BQSR {prx}.grp -o {prx}.over.bam -L {chr} -nct 4
