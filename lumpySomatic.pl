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

my ($input , $help);
GetOptions(
	"i|input=s"	=>	\$input,
	"help"	=>	\$help,
);

if ($help or ! $input){
	&help;
	exit;
}


open IN , "zcat $input | ";
my $k = 0;
while (<IN>){
	if (/^#/){
		print $_;
		next;
	}
	chomp;
	my @F = split /\t/ , $_;
	unless ($F[10] =~ /^0\/0/){
		next;
	}
	unless ($F[9] =~ /^\d\/1/){
		next;
	}
	print "$_\n";
}
close IN;

sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: lumpySomatic.pl
#
#        USAGE: ./lumpySomatic.pl  
#
#  DESCRIPTION: 
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Wang Ruiru (wangrr), wangruiru\@berrygenomics.com
# ORGANIZATION: Berry Genomics
#      VERSION: 1.0
#      CREATED: 08/05/2019 09:23:01 AM
#     REVISION: ---
#===============================================================================
EOF!
}



