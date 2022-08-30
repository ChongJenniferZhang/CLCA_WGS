echo -n "" > input.txt
for cir in */*_graph.txt
do

	base=$(basename $cir)
	dir=$(dirname $cir)
	sample=$(basename $dir)
	name=$(echo $base | sed 's/_graph.txt//')
	cys=$(echo $cir | sed 's/_graph.txt/_cycles.txt/')
	if [ -s $cys ] && [ -s $cir ]
	then
		echo -e "$name\t$cys\t$cir" >> input.txt
	fi
done

export AA_SRC=/share/work1/wangrr/local/AmpliconArchitect/src
export AA_DATA_REPO=/share/work1/wangrr/local/AmpliconArchitect/data_repo

/share/Oncology/RD_HCC/IBFC2018028/WgsReSeq/AmpliconArchitect/circle/AmpliconClassifier-main/amplicon_classifier.py --ref hg19 --input input.txt > classifier_stdout.log

