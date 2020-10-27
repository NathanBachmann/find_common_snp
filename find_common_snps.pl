#!/usr/bin/perl
#Written by Nathan Bachmann () under the supervision of Prof Vitali Sintchenko () on the 13/07/20

use warnings;
use strict;

my $usage = <<'USAGE';

USAGE:

find_common_snps.pl <trait matrix> <directory of vcf files>

trait matrix must be in tab-delimitated format

USAGE

#variable list
my $line;
my @list;
my @list2;
my $line2;
my @element;
my $file;
my %hash;
my %trait2;
my $position;
my $position2;
my @element2;
my $number;
my $number2;
my $file2;
my $line3;
my @element3;
my @element4;
my $unique1 = "0";
my $unique2 = "0";
my $cutoff1;
my $cutoff2;
my @part;
my %total1;
my %total2;
my $trait2cut;
my $traint2cut2;

#print the usage if no files are entered
unless (defined ($ARGV[0]))
{
	print $usage;
	exit;
}

#main code
my $input = $ARGV[0];
my $path = $ARGV[1];

chomp $path;

open (IN, "$input") or die "Couldn't open $input\n";

open (LOG, ">common_snp_log.txt");

#parser the input matrix trait file and store ID into 2 arrays 
while ($line = <IN>) # Need to add a check feature to make sure the trait file is formatted correctly
{
	chomp $line;
	if ($line =~ /^\w+/)
	{
		@part = split(/\t/, $line);
		if (($part[1] == "1") and ($part[2] == "0"))
		{
			push @list, $part[0];
		}
		if (($part[1] == "0") and ($part[2] == "1"))
		{
			push @list2, $part[0];
		}
		if (($part[1] == "1") and ($part[2] == "1"))
		{
			print "error: $line has 1s for both traits\n";
			exit;
		}
		if (($part[1] == "0") and ($part[2] == "0"))
		{
			print "error: $line has 0s for both traits\n";
			exit;
		}
	}
}

close IN;

#check the number of ID per list/trait
my $count = scalar(@list);
print "There are $count with trait 1\n";
print LOG "There are $count with trait 1\n";
my $count2 = scalar(@list2);
print "There are $count2 with trait 2\n";
print LOG "There are $count2 with trait 2\n";


# calculate the 80% cutoff for trait 1
$cutoff1 = $count*0.8;

# calculate the 20% cutoff for trait 2
$cutoff2 = $count2*0.2;

# calculate the 80% cutoff for trait 2 and 20% cutoff for trait 1
$trait2cut = $count2*0.8;
$traint2cut2 = $count*0.2;

print LOG "trait 1 cutoff is $cutoff1 and trait 2 cutoff is $cutoff2\n";

#print "$path\n";
print "Reading VCF files for trait 1\n";
#run these get VCF file (trait one) and stores each SNP in hash with position as the key and list of IDs
#as the value
foreach (@list)
{
	$file = $_;
	#print "$_\n";
	print LOG "reading trait 1 vcf: $file.vcf\n";
	open (IN2, "$path/$file.vcf") or die "Couldn't open $path/$file.vcf\n";
	
	while ($line2 = <IN2>)
	{ 
		chomp $line2;
		if ($line2 =~ /^\w+/)
		{
			@element = split(/\t/, $line2);
			#print "$element[3]\t$element[4]\t$file\n";
			if (($element[3] =~ /^\w{1}$/) and ($element[4] =~ /^\w{1}$/))
			{
				#print "$element[1]\t$element[3]\t$element[4]\t$file\n";
				#print "$file\t$line2\n";
				$position = $element[1];
				$total1{$position}++;  #counts the number of times this SNP appears in trait 1 set for 80/20 calculations
				if (exists $hash{$position})
				{
					$hash{$position} .= "\t$file $element[4]";
				}
				else
				{
					$hash{$position} = "$file $element[4]";
				}
			}
		}
		else
		{
			next;
		}
	}
	
	close IN2;
}

print "Reading VCF files for trait 2\n";
#extract the snp positions from the sequences in trait 2 list
foreach (@list2)
{
	$file2 = $_;
	#print "$_\n";
	print LOG "reading trait 2 vcf: $file2.vcf\n";
	open (IN3, "$path/$file2.vcf") or die "Couldn't open trait 2 file: $path/$file2.vcf\n"; 

	while ($line3 = <IN3>)
	{
		chomp $line3;
		if ($line3 =~ /^\w+/)
		{
			@element3 = split(/\t/, $line3);
			if (($element3[3] =~  /^\w{1}$/) and ($element3[4] =~ /^\w{1}$/))
			{
				$position2 = $element3[1];
				$total2{$position2}++; #counts the number of times this SNP appears in trait 2 set for 80/20 calculations
				if (exists $trait2{$position2})
				{
					$trait2{$position2} .= "\t$file2 $element3[4]";
				}
				else
				{
					$trait2{$position2} = "$file2 $element3[4]";
				}
			}
		}
		else
		{
			next;
		}
	}
	
	close IN3;
}

open (OUT, ">unique_snps_in_trait1.tab");

#print LOG "trait 1 SNPs\n";
#prints out the trait 1 hash
foreach my $key (sort keys %hash)
{
	#print "$key\t$total1{$key}\n";
	@element2 = split(/\t/,$hash{$key});
	$number = scalar(@element2); 
	#print "$key\t$number\n";
	if ($number == $count) # print out SNPs that are found in all "trait 1" genomes. 
	{
		#print "$key\t$number\t$hash{$key}\t$cutoff1\n";
		if (exists $trait2{$key})
		{
			#print "$key is not a unique SNP\n";
			next;
		}
		else
		{
			#print "$key is a true unique SNP\n";
			print OUT "$key\t$number\t$hash{$key}\n";
			$unique1++;
		}
	}
	else
	{
		next;
	}
}

close OUT;

print "Number of unique SNP with trait 1 is $unique1\n";
print LOG "Number of unique SNP with trait 1 is $unique1\n";

open (OUT2, ">unique_snps_in_trait2.tab");

#print "trait 2 SNPs\n";
foreach my $key2 (sort keys %trait2)
{
	#print "this is total check: $key2\t$total2{$key2}\n";
	@element4 = split(/\t/, $trait2{$key2});
	$number2 = scalar(@element4);
	#print "$number2\n";
	if ($number2 == $count2)
	{
		#print "$key2\t$number2\t$trait2{$key2}\n";
		if (exists $hash{$key2})
		{
			next;
		}
		else
		{
			print OUT2 "$key2\t$number2\t$trait2{$key2}\n";
			$unique2++;
		}
	}
	else
	{
		next;
	}
}

close OUT2;

print "Number of unique SNP with trait 2 is $unique2\n";
print LOG "Number of unique SNP with trait 2 is $unique2\n";

########
#code for 80/20 comparison of SNPs
########
#this section parsers the SNPs using 80/20 cutoff to be less strict
print "trait 1 cutoff is $cutoff1\n";
print "trait 2 cutoff is $cutoff2\n";

open (OUT3, ">trait1_snp_80_20.tab");

#prints out the trait 1 hash for the 80/20 comparison
foreach my $key (sort keys %hash)
{
	#print LOG "this is total match check: $key\t $total2{$key}\n";
	@element2 = split(/\t/,$hash{$key});
	$number = scalar(@element2);
	if ($number > $cutoff1) # print out trait1 SNPs that are found > cutoff "trait 1" genomes. 
	{
		if (exists $trait2{$key})
		{
			if ($total2{$key} < $cutoff2)
			{
				print OUT3 "$key\t$number\t$total2{$key}\t$hash{$key}\n";
			}
			else
			{
				next;
			}
		}
		else
		{
			print OUT3 "$key\t$number\t$hash{$key}\n";
		}
	}
	else
	{
		next;
	}
}

close OUT3;

open (OUT4, ">trait2_snp_80_20.tab");

#prints out the trait 2 hash for the 80/20 comparison
#print "trait 2 SNPs\n";
foreach my $key2 (sort keys %trait2)
{
	#print "this is total check: $key2\t$total2{$key2}\n";
	@element4 = split(/\t/, $trait2{$key2});
	$number2 = scalar(@element4);
	if ($number2 > $trait2cut)
	{
		if (exists $hash{$key2})
		{
			if ($total1{$key2} < $traint2cut2)
			{
				print OUT4 "$key2\t$number\t$total1{$key2}\t$trait2{$key2}\n";
			}
			else
			{
				next;
			}					
		}
		else
		{
			print OUT4 "$key2\t$number2\t$trait2{$key2}\n";
		}
	}
	else
	{
		next;
	}
}	

close OUT4;
close LOG;

exit;