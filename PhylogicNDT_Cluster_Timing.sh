sample=$1
mkdir -p $sample
cd $sample
/share/public/software/workspaces1/env1/PhylogicNDT/PhylogicNDT.py Cluster -i $sample  -sif sif/$sample.s1.sif -rb --maf --maf_input_type calc_ccf
/share/public/software/workspaces1/env1/PhylogicNDT/PhylogicNDT.py Timing -i $sample  -sif sif/$sample.s2.sif --driver_genes_file  driver.id
