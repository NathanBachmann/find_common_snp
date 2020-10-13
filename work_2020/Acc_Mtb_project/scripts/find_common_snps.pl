#!/usr/bin/perl
#Written by Nathan Bachmann () under the supervision of Prof Vitali Sintchenko () on the 13/07/20

use warnings;
#use strict;

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

#parser the input matrix trait file and store ID into 2 arrays 
while ($line = <IN>)
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
			print "error: $line\n";
			exit;
		}
		if (($part[1] == "0") and ($part[2] == "0"))
		{
			print "error: $line\n";
			exit;
		}
	}
}

close IN;

#check the number of ID per list/trait
my $count = scalar(@list);
print "There are $count with trait 1\n";
my $count2 = scalar(@list2);
print "There are $count2 with trait 2\n";


#print "$path\n";
print "Reading VCF files for trait 1\n";
#run these get VCF file (trait one) and stores each SNP in hash with position as the key and list of IDs
#as the value
foreach (@list)
{
	$file = $_;
	#print "$_\n";
	#print "$file.vcf\n";
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
				#print "$file\t$line2\n";
				$position = $element[1];
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
	#print "$file2.vcf\n";
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

#print "trait 1 SNPs\n";
#prints out the trait 1 hash
foreach my $key (sort keys %hash)
{
	#print "$key\t$hash{$key}\n";
	@element2 = split(/\t/,$hash{$key});
	$number = scalar(@element2); 
	#print "$number\n";
	if ($number == $count) # print out SNPs that are found in all "trait 1" genomes
	{
		#print "$key\t$number\t$hash{$key}\n";
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

open (OUT2, ">unique_snps_in_trait2.tab");

#print "trait 2 SNPs\n";
foreach my $key2 (sort keys %trait2)
{
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

exit;