for ele in enhancers gc19_pc.3utr gc19_pc.5utr lncrna.ncrna mirna.mat mirna.prom smallrna.ncrna mirna.pre lncrna.prom gc19_pc.prom	
do
	mkdir -p $ele
	Rscript ./pvalue_combination.R $ele.tsv $ele/
	Rscript adjust.R $ele/combined_p_values.automatic_method_removal.txt $ele/$ele.combined_p_values.xls
done
