#!/usr/bin/perl
BEGIN{
	push @INC,'/home/terpstramm/software/vcftools_0.1.9/lib/perl5/site_perl/';
}
use warnings;
use strict;
use Vcf;
use Data::Dumper;
my $use =<<"END";

use: perl $0 in.vcf
uses the vcftools Vcf.pm library so please add it to PERL5LIB enviroment variable
practical use example:
This dumps some basic fields:
this dumps all info/format fields + some basic fields (CHROM POS REF ALT ID QUAL FILTER)
java -jar GenomeAnalysisTK.jar \\
     -R reference.fasta \\
     -T VariantsToTable \\
     -V file.vcf \\
     -AMD \\
     -raw \\
     -F CHROM -F POS -REF -ALT -F ID -F QUAL -F FILTER \$(perl $0 file.vcf) \\
     -o results.table
for vcftools merge compatibility, i have added functionality to skip any field matching with /_[0123456789]$/
note: don't worry anymore about fixing typoos, missing any fields, exporting unpopulated fields in your VcfTableExport.

END
#

die "Invalid vcf file specified!! '$ARGV[0]'\n $use\n$!" if(not(-e $ARGV[0]));
AccumulateFieldsVCF($ARGV[0]);

sub AccumulateFieldsVCF{
	my $vcfin = $_[0];#'/home/terpstramm/Downloads/s12-193.annotated.vcf';
	my $out;
	#if($_[1] && not(-e $_[1])){
	#	my $vcfout = $_[1];
	#	open($out,'>',$vcfout) or die 'Cannot write vcf file';
	#}else{
	#	warn 'warning: no out.vcf file specified printing to STDOUT'."\n";
		$out =*STDOUT;
	#}
	
	
	my $vcf = Vcf->new(file=>$vcfin) or die "cannot open vcf file $vcfin\n";
	$vcf->parse_header();
	#print {$out} $vcf->format_header();
	$vcf->recalc_ac_an(0);
	
	my %FieldInfo;
	my $q;
	while (my $x=$vcf->next_data_hash()){
		#print Dumper($x);
		
		for my $F (keys(%{$x->{'INFO'}})){
			next if($F =~ /_[0123456789]$/);
			$FieldInfo{'F'}{$F}++;
		}
		for my $GF (@{$x->{'FORMAT'}}){
			$FieldInfo{'GF'}{$GF}++;
		}
		#print {$out} $vcf->format_line($x);
		
		
#		 if ($. > 2);
		
	}
	#
	##results here
	#
	$q = \%FieldInfo;
	#print Dumper($q);
	print " -F ".join(" -F ",sort(keys(%{$q->{'F'}})));
	print " -GF ".join(" -GF ",sort(keys(%{$q->{'GF'}})));
	$vcf->close();
}
