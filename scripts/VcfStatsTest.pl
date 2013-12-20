#!/usr/bin/perl
BEGIN{
	push @INC,'/home/terpstramm/software/vcftools_0.1.9/lib/perl5/site_perl/';
}
use warnings;
use strict;
use VcfStats;
use Vcf;
use Data::Dumper;

my @affectedSamples = ('S12_193_3');
my @parentsAffectedSamples = ('S12_193_2');
#&AddRecessiveModelInfo('/home/terpstramm/Downloads/s12-193.annotated.vcf','',\@affectedSamples,\@parentsAffectedSamples);
AccumulateFieldsVCF('/home/terpstramm/Documents/gccclusterdata/data/projects/privsmeta/forGerard/s12-204.nofilter.snps.vcf');
#RemoveFilter('/home/terpstramm/Downloads/s12-193.annotated.vcf');

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

sub AccumulateFieldsVCF{
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
	#print {$out} $vcf->format_header();
	$vcf->recalc_ac_an(0);
	
	my %FieldInfo;
	my $q;
	while (my $x=$vcf->next_data_hash()){
		#print Dumper($x);
		
		for my $F (keys(%{$x->{'INFO'}})){
			$FieldInfo{'F'}{$F}++;
		}
		for my $GF (@{$x->{'FORMAT'}}){
			$FieldInfo{'GF'}{$GF}++;
		}
		#print {$out} $vcf->format_line($x);
		
		
#		 if ($. > 2);
		
	}
	$q = \%FieldInfo;
	#print Dumper($q);
	print " -F ".join(" -F ",keys(%{$q->{'F'}}));
	print " -GF ".join(" -GF ",keys(%{$q->{'GF'}}));
	$vcf->close();
}

sub AddRecessiveModelInfo{
	#input examples			'example.vcf',	'example.recessive.out.vcf',	\@affectedSamples,	\@parentsAffectedSamples
	#or output to stdout	'example.vcf',	''							,	\@affectedSamples,	\@parentsAffectedSamples
	#minimal				'example.vcf',	''							,	\@affectedSamples
	#$x is read from a merged VCF file containing samples and controls. 
	#the controls not specified as parents are assumed not related directly.
	#parents can be omitted then the genotypes will not be required as hetrozygous.
	
	#
	##input handling
	#
	my $vcfin = $_[0];#'/home/terpstramm/Downloads/s12-193.annotated.vcf';
	my $out;
	if($_[1] && not(-e $_[1])){
		my $vcfout = $_[1];
		open($out,'>',$vcfout) or die "Cannot write vcf file:'$vcfout'\n$!";
	}else{
		warn 'warning: not out vcf file specified printing to STDOUT'."\n";
		$out =*STDOUT;
	}
	my $affectedSamples = $_[2];
	ValidateCommandlineSamplesToVcf($vcfin,$affectedSamples);

	my $parentsAffectedSamples= $_[3] if($_[3]);
	ValidateCommandlineSamplesToVcf($vcfin,$parentsAffectedSamples);
	
	my $vcf = Vcf->new(file=>$vcfin) or die "cannot open vcf file $vcfin\n";
	$vcf->parse_header();
	my (@samples) = $vcf->get_samples();
	
	my @controlSamples = @samples;
	
	SubstractArrayOfStringsFromArrayArrayOfStrings($affectedSamples, \@controlSamples);
	SubstractArrayOfStringsFromArrayArrayOfStrings($parentsAffectedSamples, \@controlSamples)if($parentsAffectedSamples);
	
	#
	##add recessive info flag header line
	#
	my $descr = "Variants found passing a recessive model for inheritance with the following groups. All samples present(".join(',',@samples)."). Defined affected samples:(".join(',',@{$affectedSamples}).")";
	if($parentsAffectedSamples){
		$descr = $descr."unaffected samples that are parents or children(".join(',',@{$parentsAffectedSamples}).")";
	}
	$vcf->add_header_line({key=>'INFO', ID=>'Recessive',Number=>0,Type=>'Flag',Description=>$descr});
	#
	##print header
	#
	print {$out} $vcf->format_header();
	
	#$vcf->recalc_ac_an(0);
	
	#
	##edit and print variant positions while 
	#
	my %is_empty_hash;
	while (my $x=$vcf->next_data_hash()){
		die if ($. > 2);
		my $allele = GetSharedHomAllele($vcf, $x, \@samples);
		TestNoHomAllele($vcf, $x, \@samples, $allele);
		TestHetsAllele($vcf, $x, \@samples, $allele);
		my ($alleles,$seps,$is_phased,$is_empty) = $vcf->parse_haplotype($x,'S12_193_2');
		#print Dumper($alleles);
		#print Dumper($seps);
		#print Dumper($is_phased);
		#print Dumper($is_empty);
		#print {$out} $vcf->format_line($x);
		if($is_empty eq 0){
			print Dumper($alleles);
			print Dumper($seps);
			print Dumper($is_phased);
			print {$out} $vcf->format_line($x);
			die if ($. > 2);
		}
		$is_empty_hash{'x'.$is_empty}++;
		
	}
	print Dumper(%is_empty_hash);
	$vcf->close();	
}

sub ValidateCommandlineSamplesToVcf{
	my $vcfin = $_[0];
	my $commandlineSamples = $_[1];
	my $vcf = Vcf->new(file=>$vcfin) or die "cannot open vcf file $vcfin\n";
	$vcf->parse_header();
	my (@samples) = $vcf->get_samples();
	$vcf->close();
	for my $commandlineSample (@{$commandlineSamples}){
		my $commandlineSampleValidated;
		for my $vcfSamples (@samples){
			$commandlineSampleValidated = 1 if($commandlineSample eq $vcfSamples);
		}
		die "invalid sample '$commandlineSample' not found in VCF samples:\n'".join("'\n'",sort(@samples))."'\n$!" if (not ($commandlineSampleValidated));
	} 
	return;
}


sub ReprintVCF{
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
		print {$out} $vcf->format_line($x);
		die if ($. > 2);
	}
	$vcf->close();
}
sub ReprintVCFnoPL{
	my $vcfin = $_[0];#'/home/terpstramm/Downloads/s12-193.annotated.vcf';
	my $out;
	if($_[1] && not(-e $_[1])){
		my $vcfout = $_[1];
		open($out,'>',$vcfout) or die 'Cannot write vcf file';
	}else{
		warn 'warning: not output vcf file specified printing to STDOUT'."\n";
		$out =*STDOUT;
	}
	
	
	my $vcf = Vcf->new(file=>$vcfin);
	$vcf->parse_header();
	
	$vcf->remove_header_line(key=>'FORMAT', ID=>'PL');
	print {$out} $vcf->format_header();
	#$vcf->recalc_ac_an(0);
	while (my $x=$vcf->next_data_hash()){
		$vcf->remove_format_field($x,'PL');
		warn Dumper($x);
		print {$out} $vcf->format_line($x);
		die if ($. > 2);
	}
	$vcf->close();
}	
sub test1 {
my $vstats = VcfStats->new(file=>'/home/terpstramm/Downloads/s12-193.annotated.vcf');
my $h =$vstats->parse_header();
#print Dumper($vstats);
#print Dumper($h);
while (my $x=$vstats->next_data_hash()){
	print join("\t",@{$x->{'FORMAT'}})."\n";
	print join("#",%{$x->{'gtypes'}})."\n";
	print join("#",keys(%{$x->{'gtypes'}}))."\n";
	#why does this fail:print join("#",@{$x}{'gtypes'}{(keys(%{$x}{'gtypes'}))}."\n";
	print keys(%{$x->{'gtypes'}})."\n";
	print Dumper($x);
	my @samples =  keys(%{$x->{'gtypes'}});
	#print join("\t",	@{$x->{'gtypes'}{ @samples }{'PL'}});
	#@{$table{$country}}
	############
	############
	print join("\t",	@{${$x}{'gtypes'}{ @samples }{'PL'}});
	#print keys(%{$x->{'gtypes'}})
	print join("\t",($x->{'CHROM'},$x->{'POS'},$x->{'REF'},@{${$x}{'ALT'}}));
	for my $sample (@samples){
		print $x->{'gtypes'}{ $sample }{'PL'}."\t";
	}
	print "\n";
	
	die "done\n" if($. > 2);
	$vstats->collect_stats($x);
}

$vstats->dump();
}
sub SubstractArrayOfStringsFromArrayArrayOfStrings{
	my $arrayObjectsToSubstract = $_[0];
	my $array = $_[1];
	my  $i1 = 0;
	while($i1 < scalar(@{$array})){
		my  $i2 = 0;
		while($i2 < scalar(@{$arrayObjectsToSubstract})){
			if(${$arrayObjectsToSubstract}[$i2] eq ${$array}[$i1]){
				splice(@{$array},$i1,1);
				$i2 = -1;
			}
			$i2++;
		}		
	$i1++;
	}

}

sub GetSharedHomAllele{
	#$vcf, $x, \@samples;
	my $vcf = $_[0];
	my $x = $_[1];
	my $samples = $_[2];
	my %alleles;
	for (@{$samples}){
		my ($alleles,$seps,$is_phased,$is_empty) = $vcf->parse_haplotype($x ,$_);
		if(not($is_empty)){
			for (@{$alleles}){
				$alleles{$_}++;
			}
		}
	}
	if(scalar(keys(%alleles))==1){
		return (keys(%alleles))[0];
	}
	return;
}
sub TestNoHomAllele{
	my $vcf = $_[0];
	my $x = $_[1];
	my $samples = $_[2];
	my $allele = $_[3];
	my $test = 1;
	for (@{$samples}){
		
		my ($alleles,$seps,$is_phased,$is_empty) = $vcf->parse_haplotype($x ,$_);
		
		if(not($is_empty)){
			my $sampletest = 1;
			for(@$alleles){
				$sampletest = 0 if(not($allele eq $_));
			}
			$test = 0 if($sampletest == 1);
		}
	}
	return $test;
}
sub TestHetAllele{
	my $vcf = $_[0];
	my $x = $_[1];
	my $samples = $_[2];
	my $allele = $_[3];
}