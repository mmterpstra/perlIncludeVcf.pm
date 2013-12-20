#!/usr/bin/perl
BEGIN{
	push @INC,'/home/terpstramm/software/vcftools_0.1.9/lib/perl5/site_perl/';
}
use warnings;
use strict;
use VcfStats;
use Vcf;
use Data::Dumper;
# input test

# use
my $use = <<"END";
	$0 in.vcf out.vcf
	resets all filtering to '.' 
END
die "no in.vcf specified on command line\n$use" if(not -e $ARGV[0]);
#main
RemoveFilter($ARGV[0],$ARGV[1]);


sub RemoveFilter{
	#add_filter($$x[6],'SnpCluster'=>1,'q10'=>0);
	my $vcfin = $_[0];#'/home/terpstramm/Downloads/s12-193.annotated.vcf';
	my $out;
	if($_[1] && not(-e $_[1])){
		my $vcfout = $_[1];
		open($out,'>',$vcfout) or die 'Cannot write vcf file';
	}else{
		warn 'warning: not out vcf file specified printing to STDOUT'."\n";
		$out =*STDOUT;
	}
	
	
	my $vcf = Vcf->new(file=>$vcfin) or die "cannot open vcf file $vcfin\n";
	$vcf->parse_header();
	print {$out} $vcf->format_header();
	#$vcf->recalc_ac_an(0);
	while (my $x=$vcf->next_data_hash()){
		#the next line actually removes the filter
		@{$x->{'FILTER'}} = ('.');
		#print Dumper($x);
		
		print {$out} $vcf->format_line($x);
		die if ($. > 2);
	}
	$vcf->close();
}