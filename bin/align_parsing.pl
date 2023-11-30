#!/usr/bin/env perl

# Parses align file, calculates read depth (RD) and genome representation (GR)
# of individual contigs; outputs contig sequences in fasta format with this
# information in the fasta header: 
# >CLXConitgY (length-average_RD-GR)
#
# RD profiles along the contig sequences can also be calculated and saved 
# in a separate file.
#
# (renamed from align_parsing_2.pl, 2010-05-10)


use Getopt::Std;
use warnings;

getopt('iop');
if ($opt_i) {
  $align_file = $opt_i;    # 
} else {
  die "Missing parameter: -i align_file\n";
}
if ($opt_o) {
  $out_file = $opt_o;    # 
} else {
  $out_file = $align_file.".info";
  print "Output file not set, using \"$out_file\"\n";
}
if ($opt_p) {
  $profile_file = $opt_p;    # RD profiles will be recorded
} else {
  $profile_file = 0;     # will be used below to switch profile output off 
  print "Parameter -p not set, RD profiles will not be saved\n";
}

open (ALIGN,$align_file) or die "Could not read from $align_file\n";
open (OUT,">$out_file") or die "Could not write to $out_file\n";
if ($profile_file) {
  open (PROF,">$profile_file") or die "Could not write to $profile_file\n";
}

while ($radek = <ALIGN>) {
  if ($radek =~/^DETAILED DISPLAY OF CONTIGS/) {    # start of alignment section
      # print $radek;
    $radek = <ALIGN>;
    while ($radek =~/\*\*\* (CL\d+) {0,1}Contig (\d+) \*\*\*/) {   # parsing individual contigs
      #  print $radek;
      $contig_id = $1."Contig".$2;
      @pokryti = ();       # array ve ktere bude pro kazdou pozici v sekvenci pokryti
      $seq = "";           # consensus sekvence contigu
#      print "$contig_id\n";
#      print $radek;
      $radek = <ALIGN>;
      while ($radek =~/\:    \.    \:/) {  # -> nasleduje dalsi blok 
        @blok = ();                        # radky z aktualniho bloku alignmentu
        do {
          $radek = <ALIGN>;
          push(@blok,$radek); 
        } while (not $radek =~/^consensus/);
        $cons_line = $radek;
        $radek = <ALIGN>;
        $radek = <ALIGN>;     # toto posledni nacteni radku slouzi v podmince tohoto i nadrazeneho cyklu while
        if (not $radek) { $radek = " "; };      # aby to nehlasilo chyby na konci souboru kde ty radky muzou chybet
      
        pop(@blok); pop(@blok);            # odstrani posledni dva radky (mezera a consensus)
      
        for ($f=10;$f<=length($cons_line);$f++) {
          $suma_pozice = 0;
          if (substr($cons_line,$f,1) =~/([A,T,G,C,N])/) {
            $seq .= $1;
            foreach $cteni (@blok) {
              if (substr($cteni,$f,1) =~/[A,T,G,C,N]/) {
                $suma_pozice++;
              }
            }
            push(@pokryti,$suma_pozice);
          }
        }
      
      }
    
      $delka_cons = @pokryti;
      $soucet = 0;
      $prumer = 0;
      foreach $suma_pozice (@pokryti) { 
        $soucet += $suma_pozice;
      }
      $prumer = sprintf("%0.1f",$soucet/$delka_cons);

      print OUT ">$contig_id ($delka_cons-$prumer-$soucet)\n";
      while ($seq_line = substr($seq,0,60,"")) {
        print OUT "$seq_line\n";
      }
      if ($profile_file) {
        print PROF ">$contig_id ($delka_cons-$prumer-$soucet)\n";
        foreach $suma_pozice (@pokryti) {
          print PROF "$suma_pozice ";
        }
        print PROF "\n";
      }
      
    }
  }
}

close ALIGN;
close OUT;
if ($profile_file) {
  close PROF;
}

