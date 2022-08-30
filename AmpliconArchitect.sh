bam=$1
prx=$2

export AA_SRC=/share/work1/wangrr/local/AmpliconArchitect/src
export AA_DATA_REPO=/share/work1/wangrr/local/AmpliconArchitect/data_repo


python /share/work1/wangrr/local/AmpliconArchitect/src/amplified_intervals.py --bed $prx.amp.bed --bam $bam --out $prx.amp.o
python /share/work1/wangrr/local/AmpliconArchitect/src/AmpliconArchitect.py --bed $prx.amp.o.bed --bam $bam --out $prx.out
