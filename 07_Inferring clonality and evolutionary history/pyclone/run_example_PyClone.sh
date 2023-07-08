
sample=sample
PyClone setup_analysis --in_files $sample.in --working_dir ./ --tumour_contents 0.11 --samples $sample --density pyclone_binomial --num_iters 1000 --prior major_copy_number --init_method connected
PyClone run_analysis --seed 5 --config_file config.yaml
PyClone build_table --config_file config.yaml --out_file $sample.pyclone_cluster.txt --table_type old_style --max_clusters 50
