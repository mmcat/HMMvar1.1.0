#/usr/bin/perl -w

use strict;

if($ARGV[0] ne "-s" ||$ARGV[2] ne "-v"||$ARGV[4] ne "-g" ||$ARGV[6] ne "-c"){
	print "Usage: perl HGVS2HMMvar.pl -s <fasta file> -v <variant file> -g <group column index> -c <HGVE allele column index>\n";
	exit;
}
# open sequence file (e.x. tp53_01_na.fas)
open(SEQ,$ARGV[1]) or die("Cannot open file for reading:$!\n");
my $seq;
my $head = <SEQ>;
while(<SEQ>){
	chomp;
	$seq.= $_;
}
# open gene variation file (e.x. my_tp53_r16_exon_only)
open(IN,$ARGV[3]) or die("Cannot open file for reading:$!\n");
my $id = 0;
my $gindex = $ARGV[5];
my $cindex = $ARGV[7];
while(<IN>){
	chomp;
	my @line = split('\t',$_);
	my ($wt,$mt,$var,$pos);
	# SNP
	if($line[$cindex]=~/c\.([0-9]+)([A-Z])>([A-Z])/){
		$pos = $1;
		$wt = $2;
		$mt = $3;
		$var = "$wt/$mt";
		print "$line[$gindex]\t$var\t$pos\n";
		$id+=1;

	}
	#del
	elsif($line[$cindex]=~/c\.([0-9]+).*del([0-9]+)/){
		$pos = $1;
		$wt = substr($seq,$1-1,$2);
		$mt = "-";
		$var = "$wt/$mt";
		print "$line[$gindex]\t$var\t$pos\n";
		$id+=1;
	}
	#ins dont know the inserted bases from the file
	elsif($line[$cindex]=~/c\.([0-9]+)_[0-9]+ins([0-9]+)/){
		$pos = $1;
		$wt = "-";
		$mt = 'X' x $2;
		$var = "$wt/$mt";
		print "$line[$gindex]\t$var\t$pos\n";
		$id+=1;
	}
	#ins know the inserted bases from the file
	elsif($line[$cindex]=~/c\.([0-9]+)_[0-9]+ins([A-Z]+)/){
		$pos = $1;
		$wt = "-";
		$mt = $2;
		$var = "$wt/$mt";
		print "$line[$gindex]\t$var\t$pos\n";
		$id+=1;
	}
	#indel
	elsif($line[$cindex]=~/c\.([0-9]+).*del([0-9]+)ins([A-Z]+)/){
		$pos = $1;
		$wt = substr($seq,$1-1,$2);
		$mt = $3;
		$var = "$wt/$mt";
		print "$line[$gindex]\t$var\t$pos\n";
		$id+=1;
	}
	else {next;}

}
