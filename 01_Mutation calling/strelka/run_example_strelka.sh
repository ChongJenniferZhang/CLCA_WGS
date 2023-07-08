{strelka} --normalBam={Nbam} --tumorBam={Tbam} --referenceFasta={genome} --runDir={dir}
./runWorkflow.py -m local -j {num}
{pipeDir}/FilterStrelka.pl -s {dir}/results/variants/somatic.snvs.vcf.gz -i {dir}/results/variants/somatic.indels.vcf.gz > {dir}/somatic.strelka.vcf 2> {dir}/som
