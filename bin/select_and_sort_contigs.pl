#!/usr/bin/env perl

# Selects contigs exceeding the specified RD 
# and sort them based on their properties
# which have to be provided in FASTA ID line:
# >CLXContigY (lenght-RD-lenght*RD)
# (RD = average Read Depth of the contig)
# (GR = genome representation (also marked as "pxd") = lenght*RD)
#
# ver. 03 - requires input filename {fn}.info.fasta
#         - output saved to {fn}.info.minRDxxx.fasta
# 
use warnings;

$jm_contigs = $ARGV[0];
if ($ARGV[1]) {
  $min_pokryti = $ARGV[1];
} else {
  die "Missing argument - minRD\n";
}
if ($jm_contigs =~/(\S+).fasta$/) {
  $jm_vystup_zakladni = $1.".minRD".$min_pokryti;
  } else {
  $jm_vystup_zakladni = $jm_contigs.".minRD".$min_pokryti;
  }
$jm_vystup_sort_RD    = $jm_vystup_zakladni."_sort-RD.fasta";       # min. read depth
$jm_vystup_sort_delka = $jm_vystup_zakladni."_sort-length.fasta";
$jm_vystup_sort_pxd   = $jm_vystup_zakladni."_sort-GR.fasta";

%pokryti = ();         # prum. pokryti contigu ctenimi
%delka_cont = ();      # delka contigu (jeho consensu)
%pxd = ();             # pxd = pokryti x delka
%sekvence = ();        # sekvence contigu

open (CONT, $jm_contigs) or die "Cannot open $jm_contigs\n";
while ($radek = <CONT>) {
  if ($radek =~/^>(CL\d+Contig\d+) \((\d+)\-(\S+)\-(\d+)\)/) {
    if ($3>=$min_pokryti) {
#      print "$1 : delka $2, prum. pokryti $3, pxd $4\n";
      $delka_cont{$1}   = $2;
      $pokryti{$1}      = $3;
      $pxd{$1}          = $4; 
    }
  }
}
close CONT;

# ze souboru contigu se nactou sekvence tech, vybranych v predchozi fazi
open (CONT, $jm_contigs) or die "Cannot open $jm_contigs\n";
$ukladat = 0;
while ($radek = <CONT>) {
  if ($radek =~/>(CL\d+Contig\d+)/) {
      if (exists($pokryti{$1})) {
          $ukladat = 1;
          $jmeno = $1;
      } else {
          $ukladat = 0;
      }
  } elsif ($ukladat) {
      $sekvence{$jmeno} .= $radek;
  }
}
close CONT;

# vygeneruji se soubory, kde jsou contigy co prosly pres >= $min_pokryti setrideny
# podle pokryti, delky, nebo delky x pokryti
open (VYSTUP, ">$jm_vystup_sort_RD") or die "Cannot write to $jm_vystup_sort_RD\n";
foreach $klic (sort {$pokryti{$b} <=> $pokryti{$a}} keys(%pokryti)) {
  print VYSTUP ">$klic ($delka_cont{$klic}\-$pokryti{$klic}\-$pxd{$klic})\n";
  print VYSTUP "$sekvence{$klic}";
}
close VYSTUP;
open (VYSTUP, ">$jm_vystup_sort_delka") or die "Cannot write to $jm_vystup_sort_delka\n";
foreach $klic (sort {$delka_cont{$b} <=> $delka_cont{$a}} keys(%delka_cont)) {
  print VYSTUP ">$klic ($delka_cont{$klic}\-$pokryti{$klic}\-$pxd{$klic})\n";
  print VYSTUP "$sekvence{$klic}";
}
close VYSTUP;
open (VYSTUP, ">$jm_vystup_sort_pxd") or die "Cannot write to $jm_vystup_sort_pxd\n";
foreach $klic (sort {$pxd{$b} <=> $pxd{$a}} keys(%pxd)) {
  print VYSTUP ">$klic ($delka_cont{$klic}\-$pokryti{$klic}\-$pxd{$klic})\n";
  print VYSTUP "$sekvence{$klic}";
}
close VYSTUP;

$celkem_pxd = 0;
foreach $hodnota (values(%pxd)) {
#  print "H: $hodnota\n";
  $celkem_pxd += $hodnota;
}
print "\nTotal GR: $celkem_pxd\n\n"; 
