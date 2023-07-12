# CLCA_WGS
This repository contains software and analysis code used for the publication: Deep Whole Genome Analysis of 494 Hepatocellular Carcinomas.

## Contents

Category | Function | Software | Version | Container name | 
--- | --- | --- | --- |--- |
1|Mutation calling|	BWA|	0.7.12|	mutation_calling.sif|
1|Mutation calling|	FASTP|	0.13.1|	mutation_calling.sif|
1|Mutation calling|	Mutect2|	4.0.11.0|	mutation_calling.sif|
1|Mutation calling|	sambamba|	0.6.8|	mutation_calling.sif|
1|Mutation calling|	samtools|	1.4|	mutation_calling.sif|
1|Mutation calling|	Strelka|	2.8.4|	strelka.sif|
2|	Identification of candidate drivers|	ActiveDriverWGS|	1.1.1|	driver.sif|
2|	Identification of candidate drivers|	dndscv|	0.1.0|	driver.sif|
2|	Identification of candidate drivers|	MutSig2CV_NC|	Only one version|	driver.sif|
2|	Identification of candidate drivers|	NBR|	Only one version|	driver.sif|
2|	Identification of candidate drivers|	pvalue_combination|	Only one version|	driver.sif|
2|	Identification of candidate drivers|	OncodriveFML|	2.3.0|	oncodrivefml.sif|
3|	Mutational signature extraction and assignment|	mSigHdp|	1.1.2|	signature.sif|
3|	Mutational signature extraction and assignment|	SigProfilerExtractor|	1.1.0|	signature.sif|
3|	Mutational signature extraction and assignment|	mSigAct|	2.2.3|	signature.sif|
4|	SV & clustered mutational events|	maftools|	2.6.05|	SV_CNA_ecDNA.sif|
4|	SV & clustered mutational events|	Shatterseek| 	0.4|	SV_CNA_ecDNA.sif|
4|	SV & clustered mutational events|	lumpy|	0.2.13|	SV_CNA_ecDNA.sif|
5|	CNA calling|	purple|	2.34|	SV_CNA_ecDNA.sif|
5|	CNA calling|	Sequenza| 	2.1.1|	SV_CNA_ecDNA.sif|
6|	Detection of ecDNA|	AmpliconClassifier| 	0.2.5|	SV_CNA_ecDNA.sif|
7|	Inferring clonality and evolutionary history|	MutationTimeR|	0.99.3|	clonality_RNA.sif|
7|	Inferring clonality and evolutionary history|	PhylogicNDT|	Only one version|	clonality_RNA.sif|
7|	Inferring clonality and evolutionary history|	PyClone| 	0.13.1|	clonality_RNA.sif|
7|	Inferring clonality and evolutionary history|	Timing_and_Signatures|	Only one version|	clonality_RNA.sif|
8|	RNA data analysis|	RSEM|	1.3.3|	clonality_RNA.sif|
8|	RNA data analysis|	STAR|	2.7.3a|	clonality_RNA.sif|










## Contact
If you have any questions or suggestions please contact us:
Chong Zhang: zhang_c@pku.edu.cn
