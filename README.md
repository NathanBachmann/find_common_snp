# find_common_snp

Important: Code is still in development and might be buggy

Description: Finds all SNPs in common with genomes that have a specific trait 
and absent in genomes without the trait. Requires a trait matrix in tab-delimited format
and a directory of vcf files

USAGE:<br>

<code>perl find_common_snps.pl <trait matrix> <directory of vcf files> [options] </code>

Input:
<trait matrix> = trait matrix must be in tab-delimitated format
<directory of vcf files> = path to a directory of vcf files

Options:
 -u {upper cutoff}: The minimum percentage of genomes in trait 1 that must have shared SNP {default: 0.8}
 -l {lower cutoff}: the maximum percentage of genomes in trait 2 that can have SNP found in trait1 {default: 0.2}

USAGE
