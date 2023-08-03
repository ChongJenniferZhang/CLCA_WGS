#!/usr/bin/env perl
use strict;
use warnings;
use utf8;
use FindBin '$Bin';
use Getopt::Long;
use Cwd 'abs_path';
use File::Basename;
use File::Path qw( mkpath );
use lib $Bin;

my ($snp , $indel , $help);
my $T_depth =8;
my $T_depth_Var =4;
my $N_depth =8;
my $N_depth_Var =2;
my $T_Var_Strand =2;
my $T_Var_Frac =0.03;
my $N_Var_Frac =0.01;
my $somaticp =0.05;

GetOptions(
	"s|snp=s"	=>	\$snp,
	"i|indel=s"	=>	\$indel,
	"help"	=>	\$help,
);

if ($help or ! $snp){
	&help;
	exit;
}

my %pon;
open IN , "zcat Wgs.PoN.vcf.gz|";
while (<IN>){
	next if /^#/;
	chomp;
	my @F = split /\t/ , $_;
	$pon{"$F[0]\t$F[1]"} = 1;
}
close IN;

my %pv = ('A'=>4,'C'=>5,'G'=>6,'T'=>7);
my $fn = 0;
for my $file ($snp , $indel){
	$fn++;
	if ($file =~ /\.gz/){
		open IN , "zcat $file|";
	}else{
		open IN , "$file";
	}
	while (<IN>){
		if (/^##/){
			print $_ if $fn == 1;
			next;
		}elsif (/^#/){
			print '#' , join("\t" , qw/CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT	TUMOR	NORMAL/) , "\n" if $fn == 1;
			next;
		}
		chomp;
		my @F = split /\t/ , $_;
		#next unless $F[6] eq 'PASS';
		next if exists $pon{"$F[0]\t$F[1]"};
		my @param=split(/;/,$F[7]);
		my $p=0;
		my $nt = '';
		my $sgt = '';
		while($p<=$#param){
			my @info=split(/[=,]/,$param[$p]);
			if ($info[0] eq 'NT'){
				$nt = $info[1];
			}elsif ($info[0] eq 'SGT'){
				$sgt = $info[1];
			}
			$p++;
		}
		next unless $nt eq 'ref';
		my ($sr , $sa) = split /\-\>/ , $sgt;
		next if $sr eq $sa;
		my ($sa1 , $sa2) = split // , $sa;
		my $gt = '';
		if ($sa1 eq $sa2){
			$gt = '1/1';
		}else{
			$gt = '0/1';
		}
		#DP:FDP:SDP:SUBDP:AU:CU:GU:TU    60:2:0:0:0,0:0,0:0,0:58,60      294:13:0:0:0,0:2,2:11,16:268,279
		#DP:DP2:TAR:TIR:TOR:DP50:FDP50:SUBDP50:BCN50	22:22:22,23:0,1:0,16:19.25:2.20:0.00:0.11	51:51:36,41:9,15:12,88:51.07:9.64:0.00:0.18
		@F[9,10] = @F[10,9];
		my @n=split(/:/,$F[10]);
		my @t=split(/:/,$F[9]);
		my ($nvar,$nref,$tvar,$tref,$ndp,$tdp);
		if ($fn==1){
			$nvar = $n[$pv{$F[4]}];
			$nref = $n[$pv{$F[3]}];
			$tvar = $t[$pv{$F[4]}];
			$tref = $t[$pv{$F[3]}];
		}else{
			$nvar = $n[3];
			$nref = $n[2];
			$tvar = $t[3];
			$tref = $t[2];				
		}
		$nvar =~ s/,\d+//;
		$nref =~ s/,\d+//;
		$tvar =~ s/,\d+//;
		$tref =~ s/,\d+//;
		$ndp = $nvar + $nref;
		$tdp = $tvar + $tref;
		next if $ndp <= 1 or $tdp <= 1;
		#if ($tvar < $nvar){
		#	@F[9,10] = @F[10,9];
		#	($tvar, $nvar, $tdp, $ndp, $tref, $nref) = ($nvar, $tvar, $ndp, $tdp, $nref, $tref);
		#}
		my $tvf = sprintf("%.4f" , $tvar/$tdp);
		my $nvf = sprintf("%.4f" , $nvar/$ndp);
		$F[8] = "GT:AD:AF:$F[8]";
		$F[9] = "$gt:$tref,$tvar:$tvf:$F[9]";
		$F[10] = "0/0:$nref,$nvar:$nvf:$F[10]";
		if($tdp>=$T_depth and $tvar>=$T_depth_Var and $ndp>=$N_depth and $nvf<$N_Var_Frac and $tvf>=$T_Var_Frac and $nvf<$tvf/5){
			print join("\t" , @F) , "\n";
		}elsif ($F[6] eq 'PASS'){
			print STDERR join("\t" , @F) , "\n";
		}
	}
	close IN;
}

sub help{
print << "EOF!";

EOF!
}



