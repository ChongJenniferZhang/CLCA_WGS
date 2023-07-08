bam=in.bam
prx=out
cnv=cnv.xls

export AA_SRC=/local/AmpliconArchitect/src
export AA_DATA_REPO=/local/AmpliconArchitect/data_repo

sed '1d' $cnv | perl -lane 'my $len = $F[2]-$F[1]+1; print "$F[0]\t$F[1]\t$F[2]\t$len\t$F[3]" if $F[3]>2.5 and $len>100000' > $prx.amp.bed

python /local/AmpliconArchitect/src/amplified_intervals.py --bed $prx.amp.bed --bam $bam --out $prx.amp.o
python /local/AmpliconArchitect/src/AmpliconArchitect.py --bed $prx.amp.o.bed --bam $bam --out $prx.out
