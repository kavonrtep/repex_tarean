# RepeatExplorer2 with TAREAN (Tandem Repeat Analyzer) #
-------------------------------------------------------------------------------
New version of RepeatExplorer with TAndem REpeat ANalyzer 

## Authors
Petr Novak, Jiri Macas, Pavel Neumann
Biology Centre CAS, Czech Republic


## How to use RepeatExplorer2

We recommend to use RepeatExplorer2 through Galaxy based web interface on our
public server at address
[https://galaxy-elixir.cerit-sc.cz](https://galaxy-elixir.cerit-sc.cz). This
server is provided as part of the [Elixir CZ project](https://www.elixir-czech.cz) and is maintained by CESNET
and CERIT-SC that are participants of this project. Server provide sufficient
computational resources for full scale RepeatExplorer2 or TAREAN analysis. This
Galaxy server also include all necessary tools for data preprocessing and QC and
some additional tools for repeat annotation.

Step by step protocols how to use RepeatExplorer on Galaxy server were published
in [Nature Protocols](http://dx.doi.org/10.1038/s41596-020-0400-y)

Free link to the full-text: [https://rdcu.be/b80Gr](https://rdcu.be/b80Gr)

Four protocols are included:

1. *De novo* repeat identification in a single species
2. Comparative repeat analysis in a set of species
3. Development of satellite DNA probes for cytogenetic experiments
4. Identification of centromeric repeats based on ChIP-seq data


If you prefer to use command line version of RepeatExplorer2 and TAREAN, read the instruction below. 

## Change log
this repository was trasfered from bitbucket (https://bitbucket.org/petrnovak/repex_tarean/) on 2023-11-30. 
[link](CHANGELOG.md)


  
## Installation ##
To use RepeatExplorer without installation, we recommend our freely
available Galaxy server at
[https://repeatexplorer-elixir.cerit-sc.cz](https://repeatexplorer-elixir.cerit-sc.cz).
This server is provided in frame of the ELIXIR-CZ project. The Galaxy
server includes also additional tools useful for data preprocessing, quality control
and genome annotation.

For installing command-line version, follow the instruction below:


Download source code using git command:

	git clone https://bitbucket.org/petrnovak/repex_tarean.git
	cd repex_tarean
    
We recommend to install dependencies using conda (conda can be installed using [miniconda](https://docs.conda.io/en/latest/miniconda.html)). The required environment can be prepared using command:

    conda env create -f environment.yml

Then activate the prepared environment using:
   
    conda activate repeatexplorer

In the `repex_tarean` direcory compile the source and prepare databases using:

	make

Support for 32-bit executables is required. If you are using Ubuntu distribution you can add 32-bit support by running:

    sudo dpkg --add-architecture i386
    sudo apt-get update
    sudo apt-get install libc6:i386 libncurses5:i386 libstdc++6:i386
   

To verify the installation, run RepeatExplorer clustering on the provided example data:

    ./seqclust -p -v tmp/clustering_output test_data/LAS_paired_10k.fas
    

### Running RepeatExplorer using Singularity:
Singularity container is available from singularity hub. Multiple versions of
repeatexplorer are available. See [singularity](https://github.com/repeatexplorer/repex_tarean) conainter definiton for
available versions. To run latest version use:
`library://repeatexplorer/default/repex_tarean:0.3.12-7a7dc9e`

Example of running RepeatExplorer2 using singularity:

get help:

    singularity exec library://repeatexplorer/default/repex_tarean:0.3.12-7a7dc9e seqclust --help


run clustering in working directory:


    mkdir working_dir && cd working_dir
    # get test data - fasta file with paired reads in interlaced format
    wget https://bitbucket.org/petrnovak/repex_tarean/raw/devel/test_data/LAS_paired_10k.fas
    singularity exec --bind ${PWD}:/data/ shub://repeatexplorer/repex_tarean seqclust -p -v /data/re_output /data/LAS_paired_10k.fas
    
    
output will be located in `working_dir/re_output`



### Running RepeatExplorer using Docker container
Singularity is prefered option but RepeatExplorer2 can be also run using Docker container((https://hub.docker.com/repository/docker/kavonrtep/repeatexplorer)) using following commands:


    # see help:
    docker run kavonrtep/repeatexplorer:2.3.8 /repex_tarean/seqclust --help
    # create working directory
    mkdir working_dir && cd working_dir
    # get test data - fasta file withj paired reads in interlaced format
    wget https://bitbucket.org/petrnovak/repex_tarean/raw/devel/test_data/LAS_paired_10k.fas
    # run clustering on test data in working directory
    docker run -v  ${PWD}:/data kavonrtep/repeatexplorer:2.3.8 /repex_tarean/seqclust -p -v /data/re_output /data/LAS_paired_10k.fas
    # output will be locate in working_dir/re_output 




## Protein databases

Repeatexplorer2 utilizes REXdb database of protein domains for repeat annotation and classification. Structure of the database is described at [http://repeatexplorer.org/](http://repeatexplorer.org/). REXdb is distributed together with RepeatExplorer source code (`database` directory)


## RepeatExplorer command line options

    usage: seqclust [-h] [-p] [-A] [-t] [-l LOGFILE] [-m {float range 0.0..100.0}]
                    [-M {0,float range 0.1..1}] [-o {float range 30.0..80.0}]
                    [-c CPU] [-s SAMPLE] [-P PREFIX_LENGTH] [-v OUTPUT_DIR]
                    [-r MAX_MEMORY] [-d DATABASE DATABASE] [-C] [-k]
                    [-a {2,3,4,5}]
                    [-tax {VIRIDIPLANTAE3.0,VIRIDIPLANTAE2.2,METAZOA2.0,METAZOA3.0}]
                    [-opt {ILLUMINA,ILLUMINA_DUST_OFF,ILLUMINA_SHORT,OXFORD_NANOPORE}]
                    [-D {BLASTX_W2,BLASTX_W3,DIAMOND}]
                    sequences
    
    RepeatExplorer:
        Repetitive sequence discovery and clasification from NGS data
    
        
    
    positional arguments:
      sequences
    
    optional arguments:
      -h, --help            show this help message and exit
      -p, --paired
      -A, --automatic_filtering
      -t, --tarean_mode     analyze only tandem reapeats without additional classification
      -l LOGFILE, --logfile LOGFILE
                            log file, logging goes to stdout if not defines
      -m {float range 0.0..100.0}, --mincl {float range 0.0..100.0}
      -M {0,float range 0.1..1}, --merge_threshold {0,float range 0.1..1}
                            threshold for mate-pair based cluster merging, default 0 - no merging
      -o {float range 30.0..80.0}, --min_lcov {float range 30.0..80.0}
                            minimal overlap coverage - relative to longer sequence length, default 55
      -c CPU, --cpu CPU     number of cpu to use, if 0 use max available
      -s SAMPLE, --sample SAMPLE
                            use only sample of input data[by default max reads is used
      -P PREFIX_LENGTH, --prefix_length PREFIX_LENGTH
                            If you wish to keep part of the sequences name,
                             enter the number of characters which should be 
                            kept (1-10) instead of zero. Use this setting if
                             you are doing comparative analysis
      -v OUTPUT_DIR, --output_dir OUTPUT_DIR
      -r MAX_MEMORY, --max_memory MAX_MEMORY
                            Maximal amount of available RAM in kB if not set
                            clustering tries to use whole available RAM
      -d DATABASE DATABASE, --database DATABASE DATABASE
                            fasta file with database for annotation and name of database
      -C, --cleanup         remove unncessary large files from working directory
      -k, --keep_names      keep sequence names, by default sequences are renamed
      -a {2,3,4,5}, --assembly_min {2,3,4,5}
                            Assembly is performed on individual clusters, by default 
                            clusters with size less then 5 are not assembled. If you 
                            want need assembly of smaller cluster set *assmbly_min* 
                            accordingly
      -tax {VIRIDIPLANTAE3.0,VIRIDIPLANTAE2.2,METAZOA2.0,METAZOA3.0}, --taxon {VIRIDIPLANTAE3.0,VIRIDIPLANTAE2.2,METAZOA2.0,METAZOA3.0}
                            Select taxon and protein database version
      -opt {ILLUMINA,ILLUMINA_DUST_OFF,ILLUMINA_SHORT,OXFORD_NANOPORE}, --options {ILLUMINA,ILLUMINA_DUST_OFF,ILLUMINA_SHORT,OXFORD_NANOPORE}
      -D {BLASTX_W2,BLASTX_W3,DIAMOND}, --domain_search {BLASTX_W2,BLASTX_W3,DIAMOND}
                            Detection of protein domains can be performed by either blastx or
                             diamond" program. options are:
                              BLASTX_W2 - blastx with word size 2 (slowest, the most sesitive)
                              BLASTX_W3 - blastx with word size 3 (default)
                              DIAMOND   - diamond program (significantly faster, less sensitive)
                            To use this option diamond program must be installed in your PATH
    


## Galaxy toolshed

RepeatExplorer can be installed on Galaxy server from toolshed repository [ https://toolshed.g2.bx.psu.edu/view/petr-novak/repeatexplorer2]( https://toolshed.g2.bx.psu.edu/view/petr-novak/repeatexplorer2).

## Reproducibility
To make clustering reproducible between runs with the 
same data, environment variable PYTHONHASHSEED must be set:

	export PYTHONHASHSEED=0
    
## Disk space requirements
Large sqlite database for temporal data is created in OS-specific temp directory- usually /tmp/ 
To use alternative location, it is necessary to specify `TEMP` environment variable.

## CPU and RAM requirements

Resources requirements can be set either from command line arguments `--max-memory` and `--cpu` or
using environment variables `TAREAN_MAX_MEM` and `TAREAN_CPU`. If not set, the pipeline will use all
available resources

## How to cite

If you use RepeatExplorer for general repeat characterization in your work please cite:

 - [Novak, P., Neumann, P., Pech, J., Steinhaisl, J., Macas, J. (2013) - RepeatExplorer: a Galaxy-based web server for genome-wide characterization of eukaryotic repetitive elements from next generation sequence read. Bioinformatics 29:792-793](http://bioinformatics.oxfordjournals.org/content/29/6/792)

or

 - [Novak, P., Neumann, P., Macas, J. (2010) - Graph-based clustering and characterization of repetitive sequences in next-generation sequencing data. BMC Bioinformatics 11 :37](http://www.biomedcentral.com/1471-2105/11/378)

If you use TAREAN for satellite detection and characterization please cite:

 - [Novak, P., Robledillo, L.A.,Koblizkova, A., Vrbova, I., Neumann, P., Macas, J. (2017) - TAREAN: a computational tool for identification and characterization of satellite DNA from unassembled short reads. Nucleic Acid Research](https://doi.org/10.1093/nar/gkx257)
 
