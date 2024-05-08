#!usr/bin/perl
use strict;
use warnings;

#Modified by GLO Sept 2013 to keep sample names and separate by tab instead of space for genotypes.

#the SNP table
my $snps = $ARGV[0];
open SNPS, $snps; #snp positions contig\tposition\tSNPs

#my $scaff = $ARGV[1];
my $scaff = "/data/home/gabriele/Spruce_Banff/references/WS77111_v2/WS77111-v2_1000plus_LGs.fa.scaffold_order";
open SCAFF, $scaff; #scaffold positions scaffold\tcontig\tstart\tend

#summary of the SNP allele frequencies
open OUT, ">$snps.contig";

my %scaff=();
#read in scaffold position
while (<SCAFF>){
	my @split=split;
	push(@{$scaff{$split[0]}}, $_);
}

#list of contig postions in scaffold
my @snp=();
my @samples;
my %samples;
while (<SNPS>){
	chomp;
	my @a1 = split(/\t/, $_);
	if ($a1[0]=~/^CHROM/){
		print OUT "CHROM\tPOS";
		foreach my $i (2..$#a1){
                        $samples{$i}=$a1[$i];
                        push(@samples,$a1[$i]);
                        print OUT "\t$a1[$i]";
                }
                print OUT "\n";
	}
       push(@snp, $_);
}

#go through each start and stop position and see if the SNP lies between them
my @split2=();
foreach my $snp(@snp){
	my @split=split ' ', $snp;
	foreach my $k (@{$scaff{$split[0]}}){
		@split2=();
		@split2=split '\t', $k;
		if ($split[1] > ($split2[2]-1) && $split[1] < ($split2[3]+1)){
			my $pos=$split[1]-$split2[2] +1;
	        	my @TAB=@split;
			shift @TAB;
        		shift @TAB;
			local $" = "\t";
			$split2[1] =~ s/>//g;  
			print OUT "$split2[1]\t$pos\t@TAB\n";
		}else{next;}
	}	
}
