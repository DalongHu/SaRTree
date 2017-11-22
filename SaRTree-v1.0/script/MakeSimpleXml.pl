#!/usr/bin/perl -w
use Bio::SeqIO;
use FindBin qw($Bin);
my $in = shift;
my $date = shift;
my $xml = shift;
unless(defined $in and defined $date and defined $xml){
	print "Usage: perl $0 <In(fasta)> <Date(list)> <Out(xml)>\n";
	exit;
}
my $fasta = Bio::SeqIO->new(-file=>"$in",-format=>'fasta');
open DATE,"$date";
my %date;
while(<DATE>){
	chomp;
	my @tmp = split /\s+/;
	$date{$tmp[0]} = sprintf "%.1f",$tmp[1];
}
close DATE;
my @tmp = sort{$a<=>$b} values %date;
my $yearRange = sprintf "%.0f",($tmp[-1]-$tmp[0])/10;
$yearRange = sprintf "%.1f",$yearRange*10;
my $yearRange10 = sprintf "%.1f",$yearRange*10;
$in =~ /([^\/]+)\.fake\.fasta/;
my $tag = $1;
my $samples = '';
foreach my $sample(keys %date){
	$samples .= '<taxon id="'.$sample.'">'."\n".'<date value="'.$date{$sample}.'" direction="forwards" units="years"/>'."\n".'</taxon>'."\n";
}
my $sequences = '';
while(my $seq = $fasta->next_seq){
	my $name = $seq->display_id;
	my $char = $seq->seq;
	$sequences .= '<sequence>'."\n".'<taxon idref="'.$name.'"/>'."\n".$char."\n".'</sequence>'."\n";
}

open TMP,"$Bin/../template/template.xml";
open OUT,">$xml";
while(<TMP>){
	s/\/\/==samples==\/\//$samples/;
	s/\/\/==sequences==\/\//$sequences/;
	s/\/\/==yearRange==\/\//$yearRange/;
	s/\/\/==yearRange10==\/\//$yearRange10/;
	s/\/\/==tag==\/\//$tag/;
	print OUT "$_";
}


