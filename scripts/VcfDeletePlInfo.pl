#!/usr/bin/perl
use strict;
use warnings;
my $use = <<"END";
	calculates filtering flags for vcf files:
		Deletes PL field per sample
	use perl $0 VCFFile.vcf
		NomalSampleID =  normal tissue other samples are assumed to be tumorSamples
		VCFFile.vcf = unannotated VCF file	
END

#my $Normal = $ARGV[0];
my $VcfFile = $ARGV[0];

open(my $in,'<',$VcfFile)or die "cannot open $VcfFile\n$!";
my $line;
	
#general stuff
my %headToIndex;
my @samples;
#Parse header / define %validChroms
my $header = 1;

while($header == 1){
	$line = <$in>;
	$line =~ s/\r$|\n$|\r\n$//g;
	if($line =~ /##/){
		#SafeTyparsing
		#warn "WARN: headerline has invalid contig/formatfield: skipped header line:'$line'\n" if(not($line =~ /^##contig=<ID=\d,length=|^##contig=<ID=|^##fileformat|^##fileDate|^##reference|^##INFO=<ID=[A-EG-Z]{2}|^##FILTER|^##FORMAT/));
		#next if(not($line =~ /^##contig=<ID=\d,length=|^##contig=<ID=|^##fileformat|^##fileDate|^##reference|^##INFO=<ID=[A-EG-Z]{2}|^##FILTER|^##FORMAT/));
		print $line."\n";
	}else{
		$header = 0;
	}
}
#Print Additional info to Header;
my $customINFO = <<"END";
END
print $customINFO;

#
##Parse VCF header / define %headToIndex
#

my $head = substr($line,1);

#print '##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">'."\n";
my @head = split("\t",$head);
my $i = 0;
map{my $tmp = $_; $headToIndex{$tmp}=$i;push(@samples,$tmp)if($i>8);$i++;}(@head);

#print header
print '#'.join("\t",(@head))."\n";

#
##do variant parsing
#

my @tabSplit;	
while($line=<$in>){

	$line =~ s/\R$//g;
	@tabSplit = split("\t",$line);	
	
	#print $headToIndex{"CHROM"}."\n";
	#die "exit on line :'$.'\n$!" if ($. == 10000);#debugging 

	my %info = infoReader($tabSplit[$headToIndex{"INFO"}]);
	#print "'".join("#",@tabSplit[@headToIndex{@samples}])."'\n";
	my $a = $tabSplit[$headToIndex{"FORMAT"}];
	my @b = @tabSplit[@headToIndex{@samples}];
	my %formatInfo = formatReader(\$tabSplit[$headToIndex{"FORMAT"}],\@b,\@samples);
	#
	###
	###########youre able to recalculate/ delete tags freely here
	###
	#
	#delete $info{"AF"};
	($tabSplit[$headToIndex{"FORMAT"}], %formatInfo) = removeFormatField("PL", \%formatInfo, \$tabSplit[$headToIndex{"FORMAT"}], \@samples);
	
	print join("\t",($tabSplit[$headToIndex{"CHROM"}],$tabSplit[$headToIndex{"POS"}],$tabSplit[$headToIndex{"ID"}],$tabSplit[$headToIndex{"REF"}],$tabSplit[$headToIndex{"ALT"}],$tabSplit[$headToIndex{"QUAL"}],$tabSplit[$headToIndex{"FILTER"}],&infoPrinter(%info),&formatPrinter(\$tabSplit[$headToIndex{"FORMAT"}],\%formatInfo)))."\n";
}

sub removeFormatField{
	warn join("\t", ($_[0],%{$_[1]},join("\t",${$_[2]}),join("\t",@{$_[3]})));
	my $removeField = $_[0];
	my %formatInfo = %{$_[1]}; #
	my $format = $$_[2]; #$tabSplit[$headToIndex{"FORMAT"}]
	my @samples = @$_[3]; #\@samples
	my @formatFields = split(":",${$_[2]});
	my @formatFieldsNew=();
	for(my $formatfield,(@formatFields)){
		push(@formatFieldsNew, $formatfield) if($formatfield ne $removeField);
	}
	$format = join(":",@formatFieldsNew);
	for my $sample (@samples){
		if($formatInfo{$sample}{$removeField}){
			delete($formatInfo{$sample}{$removeField});
		}
	}
	return($format, %formatInfo);
	
}
sub SetNormalIndex {
	my $index = 0;
	for(@{$_[1]}){
		if($_ eq ${$_[0]}){
			return $index;
		}
		$index++;
	}
	die "Name of normal sample not matched(did you make a typoo?)\n$!";
}

sub max {
    my ($max, @vars) = @_;
    for (@vars) {
        $max = $_ if $_ > $max;
    }
    return $max;
}
sub min {
    my ($min, @vars) = @_;
    for (@vars) {
        $min = $_ if $_ < $min;
    }
    return $min;
}
sub Sum {
	my $sum =0;
	for (@_){
		$sum += $_;
	}
	return $sum;
}
sub Mean {
	my $count = scalar(@_);
	my $sum = 0;
	for (@_){
		$sum += $_ if($_);
	}
	my $mean=$sum/$count;
	return $mean;
}
sub infoReader {
	my $info = shift(@_);
	my @info = split(';', $info); 
	my %info;
	map{my $tmp = $_;($info{(split('=',$tmp))[0]}=(split('=',$tmp))[1])if($tmp =~ /=/);$info{$tmp} = 'true' if(not($tmp =~ /=/))}(@info);
	return %info;
}
sub infoPrinter {
	my %info = @_ or warn "Internal error has occured at filehandle line '$.'\n$!\n";
	my @info = keys(%info);

	#GATK strigency >> don't ask dont tell..
	my @usedInfoOrder;
	my @gatkInfoOrder= ('GT','AD','DP','GQ','PL');#req order(if present): GT:AD:DP:GQ:PL
	for my $gatkInfokey (@gatkInfoOrder){
		for my $infokey (@info){
			if($infokey eq $gatkInfokey){
				push(@usedInfoOrder,$infokey);
			}
		}
	}
	for my $infokey (@info){
		my $isInGatkFormat = 0;
		for my $gatkInfokey (@gatkInfoOrder){
			if($infokey eq $gatkInfokey){
				$isInGatkFormat = 1;
			}
		}
		if($isInGatkFormat == 0){
			push(@usedInfoOrder,$infokey);
		}
	}
	my $text = "";
	for my $tmp (@usedInfoOrder){
		if($info{$tmp} eq 'true' && $text eq ""){
			$text = $tmp;
		}elsif($info{$tmp} eq 'true'){
			$text = $text.';'.$tmp;
		}elsif($info{$tmp} ne 'true' && $text eq ""){
			$text = $tmp.'='.$info{$tmp};
		}elsif($info{$tmp} ne 'true'){
			$text = $text.';'.$tmp.'='.$info{$tmp};
		}	
	}
	#map{my $tmp = $_;$text = $text.';'.$tmp$info if(not($info{$tmp} eq 'true'));$text = $text.$tmp if($info{$tmp} eq 'true')}(@info);
	return $text;
}
sub formatReader {
	my @format = split(':', ${$_[0]});
	my @sampleInfos = @{$_[1]};
	my @samples = @{$_[2]};
	my $sample_count = 0;
	my %perSampleFormatData;
	for my $sample (@samples){
		$perSampleFormatData{'sampleNames'}{$sample_count} = $sample;
		$sample_count++;
		my @formatInfo = split(":",shift(@sampleInfos));
		for my $formatfield (@format){
			$perSampleFormatData{$sample}{$formatfield} = shift(@formatInfo);
		}
	}
	$perSampleFormatData{'sampleCount'} = $sample_count;
	##while(@_){
	#	my @sampleData = split(':',shift @_);
	#	my $i=0;
	#	#map{$perSampleFormatData{$sample_count}{$format[$i]}=$_;$i++}(@sampleData);
	#	map{$formatData{$format[$i]}=$_;$i++}(@sampleData);
	#	#$sample_count++;
	##}
	#die "still sampledata left while handling file line '$.':'".join("",@_)."'\nperl error'$!'\n" if(join("",@_) ne "");
	return %perSampleFormatData;
}
sub formatPrinter {
	my @format = split(':', ${$_[0]});
	my %formatData = %{$_[1]};
	#my $sample_count = 0;
	my @sampletext = ();
	my $formatText = "";
	#GATK strigency >> don't ask dont tell..
	my @usedFormatOrder;
	my @gatkFormatOrder= ('GT','AD','DP','GQ','PL');#req order(if present): GT:AD:DP:GQ:PL
	for my $gatkFormatkey (@gatkFormatOrder){
		for my $formatkey (@format){
			if($formatkey eq $gatkFormatkey){
				push(@usedFormatOrder,$formatkey);
			}
		}
	}
	for my $formatkey (@format){
		my $isInGatkFormat = 0;
		for my $gatkFormatkey (@gatkFormatOrder){
			if($formatkey eq $gatkFormatkey){
				$isInGatkFormat = 1;
			}
		}
		if($isInGatkFormat == 0){
			push(@usedFormatOrder,$formatkey);
		}
	}
	##endof GATK stringency
	my $count = 0;
	while($count < $formatData{'sampleCount'}){
		my $perSampleTekst="";

		if($formatData{ $formatData{'sampleNames'}{$count} }{'GT'} ne './.'){
			$formatText = "";#or else nonsensical format tekst
			for my $tmp (@usedFormatOrder){
				warn "tmp:".$tmp if(not(defined($tmp)));
				warn "formatData:".$formatData{ $formatData{'sampleNames'}{$count} }{$tmp} if(not(defined($formatData{ $formatData{'sampleNames'}{$count} }{$tmp})));
				if($formatData{ $formatData{'sampleNames'}{$count} }{$tmp} eq 'true' && $formatText eq ""){
					$formatText = $tmp;
					$perSampleTekst = $tmp;
				}elsif($formatData{ $formatData{'sampleNames'}{$count} }{$tmp} eq 'true'){
					$formatText = $formatText.':'.$tmp;
					$perSampleTekst = $perSampleTekst.':'.$tmp;
				}elsif($formatData{ $formatData{'sampleNames'}{$count} }{$tmp} ne 'true' && $formatText eq ""){
					$formatText = $tmp;
					$perSampleTekst = $formatData{ $formatData{'sampleNames'}{$count} }{$tmp};
				}elsif($formatData{ $formatData{'sampleNames'}{$count} }{$tmp} ne 'true'){
					$formatText = $formatText.':'.$tmp;
					$perSampleTekst = $perSampleTekst.':'.$formatData{ $formatData{'sampleNames'}{$count} }{$tmp};
				}
			}
		}else{
			$perSampleTekst='./.';
		}
		push(@sampletext,$perSampleTekst);
		$count++;
	}
	
	return ($formatText, join("\t",@sampletext));
}
#sub arrayFraction {
#	my $DP = shift(@_);
#	my @division;
#	for my $EC (@_){
#		push(@division,$EC/$DP);
#	}
#	return @division;
#}
sub calcADvals {
	my $DP = shift(@_);
	my @EC = split(",",shift(@_));
	my @ADvals;
	$ADvals[0] = $DP;
	for my $EC (@EC){
		push(@ADvals,$EC);
		$ADvals[0] =$ADvals[0] - $EC;
	}
	return @ADvals;
}