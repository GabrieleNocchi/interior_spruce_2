#!usr/bin/perl
use strict;
use warnings;

#read in fasta file
my $FASTA = "$ARGV[0]";
my $OUTNAME = "$ARGV[1]";
my $no=100000; # number of scaffold - change depending on the number you would like
my $i;

open (FASTA, "<$FASTA") || die "$_ not found.\n";

my @header=();
my @fasta=();
my @temp=();
my $a=0;

#read in fasta file into an array
while (<FASTA>) {
        chomp $_ ;
        if(/\>/){
        	push(@header, $_);
		if($a==1){
			my $join= join '', @temp;
			push (@fasta, $join);
			@temp=();
		}
		$a=1;
	}else{
		push (@temp, $_);

	}

}
my $join= join '', @temp;
push (@fasta, $join);

close FASTA;

#calculate the number of contigs/scaffold
my $size = scalar(@header);
my $F=$size/$no;
print "$size of array $size, F $F\n";

#the buffer of As between each contig
my $A='AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'; #30 A


open OUT1, ">$FASTA.pseudo";
open OUT2, ">$OUTNAME";

my $end=0;
my $c=0;

#go through each scaffold
for ($i=1; $i<=$no; $i++){
	print OUT1 ">Scaffold_$i\n";
	$end=int($F+$c);
	my @scaffolds=();
	#put the remaining contigs in the last scaffold
	if ($i==$no){
		$end=$size;
	}

	my $start=1;
	my $d=$c;

	#once large enough start a new scaffold
	for (my $j = $d; $j < $end; $j++){
		#print out contig
		push(@scaffolds, $fasta[$j]);
		my $length=length($fasta[$j]);
		my $stop=$start+$length-1;
		push (@scaffolds, $A);
		#print out scaffold, contig, start and end position in contig + length
		print OUT2 "Scaffold_$i\t$header[$j]\t$start\t$stop\t$length\n";
		$start=$stop+length($A)+1;
		$c++;
	}

	#print out 70 bases each line for the scaffolds so it is easier to read
	my $join = join "", @scaffolds;
	my @arrayBases = split ( //, $join);
	my $n=0;
	foreach( @arrayBases ) {
		print OUT1 $_;
		if( ($n + 1) % 70  == 0 ) {					# Adds a new line after every 70 bases for more readable and managable output
			print OUT1 "\n";

		}

	$n++;
	}

	print OUT1 "\n";
}

# to Add
sed '/^$/d' WS77111-v2_1000plus_LGs.pseudo.fa > temp66.txt
mv temp66.txt WS77111-v2_1000plus_LGs.pseudo.fa
