j=Chimeric.out.junction
prx=$1
n=Log.final.out
/share/public/software/Onc_Soft/STAR-Fusion-STAR-Fusion-v1.8.1/STAR-Fusion --chimeric_junction $j --genome_lib_dir fusion/lib/ctat_genome_lib_build_dir/ --output_dir $prx --CPU 8 --totalnumberreads $(cat $n | grep 'Number of input reads' | perl -lane 'print $F[-1]') --examine_coding_effect
