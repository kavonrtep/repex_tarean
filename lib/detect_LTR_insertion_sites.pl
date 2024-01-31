#!/usr/bin/env perl

# parses ACE files of assembled repeats and detects potential
# LTR borders/insertion sites of LTR-retroelements

# "site" is a region (size of $window) including TG or CA
# "out" is a region adjacent to the site, presumably representing insertion sites

# this is RepeatExplorer version of "detect_insertion_sites_LTRs.pl"
# -m default set to 10 

use Getopt::Std;




getopt('iowsmdrp');
if ($opt_i) {
	$infile = $opt_i;
} else {
	die "-i input_file_name missing\n";
}

if ($opt_p) {
    $db_PBS = $opt_p;
} else {
    die "-p PBS database is missing\n";
}




if ($opt_o) {
	$outfile = $opt_o;
} else {
	die "-o output_file_name missing\n";
}
if ($opt_w) {
	$window = $opt_w;
} else {
	$window = 7;
	print "window size not set, using default ($window)\n";
}
if ($opt_s) {
	$min_site_depth = $opt_s;   # minimal average read depth (over $window) required for the site
} else {
	$min_site_depth = 10;
	print "min_site_depth not set, using default ($min_site_depth)\n";
}
if ($opt_m) {
	$min_out_masked = $opt_m;  # minimal average number of masked reads outside the site (over $window)
} else {
	$min_out_masked = 10;
	print "min_out_masked not set, using default ($min_out_masked)\n";
}
if ($opt_d) {
	$min_masked_fold_diff = $opt_d;  # how many times should the proportion of masked reads "out" be higher than in "site"
} else {
	$min_masked_fold_diff = 3;
	print "min_masked_fold_diff not set, using default ($min_masked_fold_diff)\n";
}
if ($opt_x) {
	$max_char_to_masked = $opt_x;  # max fold difference between depth in "site" and masked depth "out"
} else {
	$max_char_to_masked = 10;
	print "max_char_to_masked not set, using default ($max_char_to_masked)\n"; 
}
if ($opt_r) {
	$extract_region = $opt_r;
} else {
	$extract_region = 30;
	print "extract_region not set, using default ($extract_region)\n";
}

# main
$out_table = $outfile;
$out_LTR   = "$outfile.LTR";
$out_ADJ   = "$outfile.ADJ";
open (IN,$infile) or die;
open (OUT,">$out_table") or die;
open (LTR,">$out_LTR") or die;  # LTR end regions as fasta seq; all are converetd to ....CA (so TG... regions are reverse-complemented)
open (ADJ,">$out_ADJ") or die;  # regions adjacent to LTR ends; if LTR end is rev-complemented, so is its corresponding adjacent region
print OUT "#Parameters:\n";
print OUT "#infile\t$infile\n#outfile\t$outfile\n#window\t$window\n#min_site_depth\t$min_site_depth\n";
print OUT "#min_out_masked\t$min_out_masked\n#min_masked_fold_diff\t$min_masked_fold_diff\n#max_char_to_masked\t$max_char_to_masked\n#extract_region\t$extract_region\n\n";
print OUT "CL\tcontig\tTG/CA\tposition\tsite\tsite_depth\tout_masked\tmasked_ratio_site\tmasked_ratio_out\tregion_in\tregion_out\tblast PBS\n";
print "Analyzing ACE file...\n";
$prev = 0;
while ($radek = <IN>) {
	$contig_found = &read_contig;
	if ($contig_found) {
		if ($cl > $prev) {
			$prev = $cl;
		}
		&reconstruct_assembly;
		&find_sites;
	}
}
close IN;
close OUT;
close LTR;
close ADJ;
print "Running blast against tRNA database...\n";
&add_PBS_info;    # detects similarities of sequences in ADJ to tRNA database (!!! reads ADJ and $out_table !!!)

$error = system("rm $out_table");
if ($error) {
	print "Error removing $out_table\n";
}

sub read_contig {
	my ($reads_found,$read_id);
	# global variables
	$cl = 0;
	$contig = 0;
	$cont_length = 0;
	$reads = 0;   # number of reads
	$cons = "";   # contig consensus (including gaps *)
	%read_starts = ();  # starts of reads within assembly
	%read_lengths = (); # length of reads in assembly (may contain gaps)
	%read_from = ();    # start of non-masked part of read sequence (relative to the read)
	%read_to = ();      # end of non-masked part of read sequence   
	
	do {
		if ($radek =~/^CO CL(\d+)Contig(\d+) (\d+) (\d+)/) {
			$cl = $1; $contig = $2; $cont_length = $3; $reads = $4;
			while ($radek = <IN> and length($radek) > 1) {
				chomp($radek);
				$cons .= $radek;
			}
			do {
				if ($radek =~/^AF (\S+) [UC] ([-]?\d+)/) {
					#print "$1 : $2\n";
					$read_starts{$1} = $2;
				}
			} while ($radek = <IN> and not $radek =~/^BS \d+ \d+/);
			$reads_found = 0;
			while ($reads_found < $reads) {        # expects previously specified number of reads 
				$radek = <IN>;              # expects RD lines with read ids alternate with QA lines
				if ($radek =~/^RD (\S+) (\d+)/) {
					$read_id = $1;
					$read_lengths{$read_id} = $2;
				}
				if ($radek =~/^QA (\d+) (\d+)/) {
					$read_from{$read_id} = $1;
					$read_to{$read_id} = $2;
					$reads_found++;
				}
			}
			return 1;
		}
	} while ($radek = <IN>);
	return 0;
}

sub reconstruct_assembly {
	my ($id,$min_start,$max_end,$shift,$f,$poz);
	# global variables
	@assembly_seq = ();     # sequence at each position of assembly; it corresponds to consensus, or regions
	                        # before (-) and after (+) consensus [assembly may be longer than consensus]
	@assembly_char = ();    # number of assembly positions with non-masked characters (nucleotides or gaps)
	@assembly_masked = ();  # number of masked positions
	$assembly_length = 0;   # total length, including - and + regions
	
	$min_start = 1; $max_end = 1;
	foreach $id (keys(%read_starts)) {
		if ($min_start > $read_starts{$id}) {
			$min_start = $read_starts{$id};
		}
		if ($max_end < $read_starts{$id}+$read_lengths{$id}-1) {
			$max_end = $read_starts{$id}+$read_lengths{$id}-1;
		}
	}
	if ($min_start < 1) {
		$shift = abs($min_start) +1;
		$assembly_length = $shift + $max_end;
	} else {
		$assembly_length = $max_end;
		$shift = 0;
	}
	
	for ($f=1;$f<=$shift;$f++) {
		$assembly_seq[$f] = "N";
	}
	for ($f=1;$f<=$cont_length;$f++) {
		$assembly_seq[$f+$shift] = substr($cons,$f-1,1);
	}
	for ($f=$shift+$cont_length+1;$f<=$assembly_length;$f++) {
		$assembly_seq[$f] = "N";
	}

	for ($f=1;$f<=$assembly_length;$f++) {
		$assembly_char[$f] = 0;
		$assembly_masked[$f] = 0;
	}
	foreach $id (keys(%read_starts)) {
		for ($f=1;$f<=$read_lengths{$id};$f++) {
			$poz = $read_starts{$id} + $shift + $f - 1;
			if ($f>=$read_from{$id} and $f<=$read_to{$id}) {
				$assembly_char[$poz]++;
			} else {
				$assembly_masked[$poz]++;
			}
		}
	}
}

sub revcompl {
	my ($input) = @_;
	my ($base,$seq,$f);
	
	$seq = "";
	for ($f=length($input)-1;$f>=0;$f--) {
		$base = substr($input,$f,1);
		if ($base eq "A") {
			$seq .= "T";
		} elsif ($base eq "T") {
			$seq .= "A";
		} elsif ($base eq "G") {
			$seq .= "C";
		} elsif ($base eq "C") {
			$seq .= "G";
		} elsif ($base eq "+" or $base eq "-") {
			$seq .= $base;
		} else {
			$seq .= "N";
		}
	}
	return $seq;
}

sub find_sites {
	my ($f,$site_sum_char,$site_sum_masked,$out_sum_char,$site_seq,$out_sum_masked,@TG,@CA,$pos);
	my ($masked_ratio_site,$masked_ratio_out,$site_depth,$out_masked,$region);
	
	# find positions of LTR borders (TG and CA)
	@TG = ();   # positions of "T"
	@CA = ();   # positions of "A"
	for ($f=1;$f<$assembly_length;$f++) {
		if ($assembly_seq[$f] eq "T" and $assembly_seq[$f+1] eq "G") {
			push(@TG,$f);
		}
		if ($assembly_seq[$f] eq "C" and $assembly_seq[$f+1] eq "A") {
			push(@CA,$f+1);
		}
	}
	
	foreach $pos (@TG) {
		if ($pos-$window > 0 and $pos+$window-1 <= $assembly_length) {
			$site_sum_char = 0; $site_sum_masked = 0; $site_seq = "";
			for ($f=$pos;$f<$pos+$window;$f++) {
				$site_sum_char += $assembly_char[$f];
				$site_sum_masked += $assembly_masked[$f];
				$site_seq .= $assembly_seq[$f];
			}
			$out_sum_char = 0; $out_sum_masked = 0;
			for ($f=$pos-$window;$f<$pos;$f++) {
				$out_sum_char += $assembly_char[$f];
				$out_sum_masked += $assembly_masked[$f];
			}
			$site_depth = sprintf("%0.1f",$site_sum_char/$window);   # average read (unmasked) depth over the site
			$out_masked = sprintf("%0.1f",$out_sum_masked/$window);  # average number of masked reads outside the site
			$masked_ratio_site = sprintf("%0.4f",$site_sum_masked/($site_sum_masked+$site_sum_char));
			$masked_ratio_out  = sprintf("%0.4f",$out_sum_masked/($out_sum_masked+$out_sum_char));
			if ($site_depth >= $min_site_depth and $out_masked >= $min_out_masked) {
				if ($masked_ratio_out >= ($min_masked_fold_diff * $masked_ratio_site) and $max_char_to_masked >= ($site_depth/$out_masked)) {
					print OUT "$cl\t$contig\tTG\t$pos\t$site_seq\t$site_depth\t$out_masked\t$masked_ratio_site\t$masked_ratio_out\t";
					$region = "";
					for ($f=$pos;$f<=$assembly_length;$f++) {
						if ($assembly_seq[$f] ne "*") {
							$region .= $assembly_seq[$f];
						}
						if (length($region) == $extract_region) {
							$f = $assembly_length;  # terminate cycle
						}
					}
					print OUT "$region\t";
					print LTR ">CL",$cl,"c".$contig."_TG_$pos\n";
					$region = &revcompl($region);
					print LTR "$region\n";
					$region = "";
					for ($f=$pos-1;$f>0;$f=$f-1) {
						if ($assembly_seq[$f] ne "*") {
							$region = $assembly_seq[$f].$region;
						}
						if (length($region) == $extract_region) {
							$f = 0;  # terminate cycle
						}
					}
					print OUT "$region\n";
					print ADJ ">CL",$cl,"c".$contig."_TG_$pos\n";
					$region = &revcompl($region);
					print ADJ "$region\n";
				}
			}
		}
	}
	
	foreach $pos (@CA) {
		if ($pos-$window+1 > 0 and $pos+$window <= $assembly_length) {
			$site_sum_char = 0; $site_sum_masked = 0; $site_seq = "";
			for ($f=$pos-$window+1;$f<=$pos;$f++) {
				$site_sum_char += $assembly_char[$f];
				$site_sum_masked += $assembly_masked[$f];
				$site_seq .= $assembly_seq[$f];
			}
			$out_sum_char = 0; $out_sum_masked = 0;
			for ($f=$pos+1;$f<=$pos+$window;$f++) {
				$out_sum_char += $assembly_char[$f];
				$out_sum_masked += $assembly_masked[$f];
			}
			$site_depth = sprintf("%0.1f",$site_sum_char/$window);   # average read (unmasked) depth over the site
			$out_masked = sprintf("%0.1f",$out_sum_masked/$window);  # average number of masked reads outside the site
			$masked_ratio_site = sprintf("%0.4f",$site_sum_masked/($site_sum_masked+$site_sum_char));
			$masked_ratio_out  = sprintf("%0.4f",$out_sum_masked/($out_sum_masked+$out_sum_char));
			if ($site_depth >= $min_site_depth and $out_masked >= $min_out_masked) {
				if ($masked_ratio_out >= ($min_masked_fold_diff * $masked_ratio_site) and $max_char_to_masked >= ($site_depth/$out_masked)) {
					print OUT "$cl\t$contig\tCA\t$pos\t$site_seq\t$site_depth\t$out_masked\t$masked_ratio_site\t$masked_ratio_out\t";
					$region = "";
					for ($f=$pos;$f>0;$f=$f-1) {
						if ($assembly_seq[$f] ne "*") {
							$region = $assembly_seq[$f].$region;
						}
						if (length($region) == $extract_region) {
							$f = 0;  # terminate cycle
						}
					}
					print OUT "$region\t";
					print LTR ">CL",$cl,"c".$contig."_CA_$pos\n";
					print LTR "$region\n";
					$region = "";
					for ($f=$pos+1;$f<=$assembly_length;$f++) {
						if ($assembly_seq[$f] ne "*") {
							$region .= $assembly_seq[$f];
						}
						if (length($region) == $extract_region) {
							$f = $assembly_length;  # terminate cycle
						}
					}
					print OUT "$region\n";
					print ADJ ">CL",$cl,"c".$contig."_CA_$pos\n";
					print ADJ "$region\n";
				}
			}
		}
	}
}

sub add_PBS_info {
	my ($pbs_blast_command,@pol,$rad,$prev_query,@table,$tab_length);
	
	$pbs_blast_command = "blastall -p blastn -d $db_PBS -i $out_ADJ -m 8 -b 1 -e 1 -W 7 -F F";
	
	@table = ();
	open (TAB,$out_table) or die;
	while ($rad = <TAB>) {
		push(@table,$rad);
		$tab_length++;
	}
	close TAB;
	
	open (BLAST,"$pbs_blast_command |") or die;
	$prev_query = "";
	while ($rad = <BLAST>) {
		if ($rad =~/^CL(\d+)c(\d+)_(TG|CA)_(\d+)\t\S+\t\S+\t/) {   
			if ("$1\t$2\t$3\t$4" ne $prev_query) {           # to exclude additional HSPs from the same query/subject pair
				for ($f=0;$f<$tab_length;$f++) {       
					@pol = split(/\t/,$table[$f]);
					if ($pol[0] eq "$1" and $pol[1] eq "$2" and $pol[2] eq "$3" and $pol[3] eq "$4") {
						chomp($table[$f]);
						$table[$f] .= "\t$rad";
						$f = $tab_length;  # terminate cycle
					}
				}
				$prev_query = "$1\t$2\t$3\t$4";
			}
		}
	}
	close BLAST;
	
	open (TAB_WITH_BLAST,">$out_table.with_PBS_blast.csv") or die;
	for ($f=0;$f<$tab_length;$f++) {
		print TAB_WITH_BLAST $table[$f];
	}
	close TAB_WITH_BLAST;
}




