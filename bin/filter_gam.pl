#!/bin/env perl

$infile = $ARGV[0];
$bamstat  = $ARGV[1];
$datadir = $ARGV[3];

$stat_file = $datadir."/probe_level_stat";
$white_list = $datadir."/white_list";
$head_file = $datadir."/data_head";

open WLIST, $white_list;
@white_list = ();
map{
	chomp;
	push @white_list, $_;
}<WLIST>;

open BAMSTAT, $bamstat;
$ctreads, $controlreads = 0, 0;
map{
	chomp;
	@line = split/\t/;
	$ctreads += $line[2] if ($line[0] >= 141 && $line[0] <= 152);
	$controlreads += $line[2] if ($line[0] >= 165 && $line[0] <= 250);
}<BAMSTAT>;

$x = $ctreads/$controlreads;
$tumorfraction = $x < 0.21 ? 0.0123 : 1.23 * ($x - 0.2);
$tumorfraction = $tumorfraction <= 0.9 ? $tumorfraction : 0.9;

$threshold = 3;
open STAT, $stat_file;
%mean = ();
%stdev = ();
map{
	chomp;
	@line = split/\t/;
	$iden = join("\t",@line[0..3]);
	$mean{$iden} = $line[4];
	$stdev{$iden} = $line[5];
}<STAT>;

open INFILE, $infile;
%count = ();
%ampprobes = ();
%zscore = ();
%sumlog2 = ();

map{
	chomp;
	@line = split/\t/;
	$iden = join("\t",@line[0..3]);
	$count{$line[3]}++;
	$zscore{$line[3]} += ($line[12] - $mean{$iden}) / $stdev{$iden};
	$sumlog2{$line[3]} += $line[12] - $mean{$iden};
	if (($line[12] - $mean{$iden}) / $stdev{$iden} > $threshold + 1) {
		$ampprobes{$line[3]} += 1;
	}
}<INFILE>;

open OUTP, $head_file;
print "chr\tstart\tend\tgene_symbol\taccession_number\tnumber_of_probes\tgt_probe_number\tindex\tthreshold_index\tcopy_number\n";
while(<OUTP>) {
	chomp;
	@line=split/\t/;
	$thr = $threshold * sqrt($count{$line[3]});
	next unless $line[3] ~~ @white_list;
	if ($zscore{$line[3]} / sqrt($count{$line[3]}) > $thr && $ampprobes{$line[3]} > $line[5] * 0.4) {
		$meanlog2 = $sumlog2{$line[3]} / $count{$line[3]};
		$copynumber = (2 ** (1 + $meanlog2) - 2 * (1 - $tumorfraction)) / $tumorfraction;
		next unless $copynumber >= 8 ;
		$copynumber = $copynumber > 20 ? 20 : int($copynumber + 0.5);
		print ($_, "\t", $ampprobes{$line[3]}, "\t", $zscore{$line[3]} / sqrt($count{$line[3]}), "\t", $thr, "\t", $copynumber, "\n");
	}
}
