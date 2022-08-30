
sample=$1
PyClone setup_analysis --in_files $sample.in --working_dir $dir/ --tumour_contents $(grep $sample purity.xls | cut -f 2) --samples $sample --num_iters 100 --prior major_copy_number --init_method connected
PyClone run_analysis --seed 5 --config_file config.yaml
PyClone build_table --config_file config.yaml --out_file $sample.pyclone_cluster.txt --table_type old_style
