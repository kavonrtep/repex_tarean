# RepeatExplorer2 Changelog

## 0.3.8.1 (2.3.8.1 version on Galaxy)
May 27, 2022

   - additional output in comparative analysis - observer/expected edges ratio is now available as CSV file

## 0.3.8. (2.3.8. version on Galaxy)
Apr 9 2020

   - Some additional output added  - links to fasta and csv files in HTML report
   - Sensistive all-to-all option added

## 0.3.7 (2.3.7 version on Galaxy)
Jan 9 2020

   - script for preparation of toolshed package added
   - RepeatExplorer2 and TAREAN are available also as Galaxy toolshed package: https://toolshed.g2.bx.psu.edu/view/petr-novak/repeatexplorer2
   - better version reporting 

##  0.3.6
Dec 20 2019
  
   - some dependencies removed 
   - all dependencies can be installed using miniconda

## v0.3.5
Nov 22 2019

   - better reporting when filtering of abundant satellite is used
   - updated format of HTML report
   - number of changes in names of output files - file names are now more informative. Changes include:
   
| old file name                                                            | new file name                                                        |
|--------------------------------------------------------------------------|----------------------------------------------------------------------|
| TR_consensus_rank_1_.fasta                                               | TAREAN_consensus_rank_1.fasta                                        |
| TR_consensus_rank_2_.fasta                                               | TAREAN_consensus_rank_2.fasta                                        |
| TR_consensus_rank_3_.fasta                                               | TAREAN_consensus_rank_3.fasta                                        |
| TR_consensus_rank_4_.fasta                                               | TAREAN_consensus_rank_4.fasta                                        |
| sequences/sequences.fasta                                                | reads/reads.fasta                                                    |
| clustering/clusters/dir_CL0001/reads.fas.CL1.aln.info.minRD5_sort-GR     | clustering/clusters/dir_CL0001/contigs.info.minRD5_sort-GR.fasta     |
| clustering/clusters/dir_CL0001/reads.fas.CL1.aln.info.minRD5_sort-length | clustering/clusters/dir_CL0001/contigs.info.minRD5_sort-length.fasta |
| clustering/clusters/dir_CL0001/reads.fas.CL1.aln.info.minRD5             | clustering/clusters/dir_CL0001/contigs.info.minRD5_sort-RD.fasta     |
| clustering/clusters/dir_CL0001/reads.fas.CL1.aln.profile                 | clustering/clusters/dir_CL0001/contigs.profile                       |
| clustering/clusters/dir_CL0001/reads.fas.CL1.aln.info                    | clustering/clusters/dir_CL0001/contigs.info.fasta                    |
| clustering/clusters/dir_CL0001/reads.fas.CL1.ace                         | clustering/clusters/dir_CL0001/contigs.ace                           |
| clustering/clusters/dir_CL0001/reads.fas.CL1.contigs.qual                | clustering/clusters/dir_CL0001/contigs.qual                          |
| clustering/clusters/dir_CL0001/reads.fas.CL1.aln                         | clustering/clusters/dir_CL0001/contigs.aln                           |
| clustering/clusters/dir_CL0001/reads.fas.CL1.contigs                     | clustering/clusters/dir_CL0001/contigs.fasta                         |
| clustering/clusters/dir_CL0001/reads.fas.CL1.info                        | clustering/clusters/dir_CL0001/assembly.info                         |
| clustering/clusters/dir_CL0001/reads.fas.CL1.singlets                    | clustering/clusters/dir_CL0001/singlets.fasta                        |
| clustering/clusters/dir_CL0001/reads.fas.CL1.contigs.links               | clustering/clusters/dir_CL0001/contigs.links                         |
| clustering/clusters/dir_CL0001/reads_oriented.fas                        | clustering/clusters/dir_CL0001/reads_selection_oriented.fasta        |
| clustering/clusters/dir_CL0001/CL1_directed_graph.RData.                 | clustering/clusters/dir_CL0001/graph_layout_directed.RData           |
| clustering/clusters/dir_CL0001/CL1_tmb.png                               | clustering/clusters/dir_CL0001/graph_layout_tmb.png                  |
| clustering/clusters/dir_CL0001/CL1.png                                   | clustering/clusters/dir_CL0001/graph_layout.png                      |
| clustering/clusters/dir_CL0001/CL1.GL                                    | clustering/clusters/dir_CL0001/graph_layout.GL                       |
| clustering/clusters/dir_CL0001/blast.csv                                 | clustering/clusters/dir_CL0001/hitsort_part.csv                      |

## v0.3.4 (1.0.0 version on Galaxy)
Oct 31 2019

  - Classification of superclusters improved, classification now uses information about LTR/PBS for classification of Class_I elements
  

## v0.3.2
Oct 9 2019

  - Graphical reporting of comparative analysis added 

## v0.3.1
Jan 9 2019

  - Improved detection of low complexity repeats and satellites with shorter monomer

## v0.3.0
Oct 25 2018

 - For back-compatibility, it is possible to select protein database version
 - Databases of protein domains went public 

## v0.2.10
Oct 24 2018

 - Protein database for Viridiplantae updated

## v0.2.9
Jan 3 2018

 - read depth for contigs is calculated
 - contigs are also sorted based on genome representation
 
## v0.2.8
Dec 18 2017

 - by default assembly is done on clusters with size at least 5 reads
 - preset option of Illumina, Illumina short and Oxford Nanopore
 - option for analyzing Metazoa/ Viridiplantae added
 - protein databases can be obtained **using fetch_databases.sh** script (password required)

## v0.2.7
Dec 05 2017

 - improved DNA database - protein domain were masked not to interfere with protein domain database
 - another bug fix in parallelization and assembly 
 - alternative to blastx added - DIAMOND program can be used instead
 
## v0.2.6
Nov 14 2017

- Improved classification of superclusters
- tables in HTML reports improved
- assembly is performed on low confidence satellite sequences
- assembly optimized, bug fix in parallelization

## v0.2.5
Sep 1 2017

+ More options added to _galaxy interface_ of TAREAN and RepeatExplorer2:
    - automatic satellite filtering
    - cluster size threshold setting
    - keep original sequence names option
+ Pipeline tests added

## v0.2.4
Aug 10 2017

- Filtering of abundant satellite sequences
- Improved (more sensitive) search of protein domains 
- bug fix in parallelization
- SHORT_ILLUMINA option added, this enable to analyze shorter reads (50nt) - command line only
- OXFORD_NANOPORE option added (experimental feature, command line only)

## v0.2.3
Jun 9 2017

- some changes in standalone TAREAN

## v0.2.2

- improved visualization of comparative analysis 

## v0.2.1

- log file is automatically copied to output and is part of archive(in galaxy)

## v0.2.0

+ Improved HTML output
+ Updated documentation
+ Full clustering analysis - RepeatExplorer2
+ New features:
    - Cap3 assembly added
    - Comparative analysis
    - LTR detection from contigs
    - Custom database option
    - Bugfix in excessive CPU usage
    - Detection of TE protein domains

## v0.1.1

+ New features:
    - Cluster merging
    - Scan of clusters against DNA database (detection of rDNA, plastid, mitochondrial and contaminants)

## v0.1.0



