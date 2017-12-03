# Geneparser
Gene parsing tool for finding shared genes across genomes  

## Prerequisites
Python 3.6+  
NCBI blastp (available in $PATH)  
NCBI makeblastdb (available in $PATH)  

## Install
- Clone this repository
    `git clone https://github.com/mmmckay/geneparser.git`
- Enter the directory
    `cd geneparser`
- Install the required Python packages
    `pip3.6 install -r requirements.txt`

## Input
Geneparser takes in 2 files per genome  
1. A rawdata .csv file containing a single genome BLASTed against all other desired genome  
    Rawdata columns are 'query', 'subject', 'Percent Identity', 'Percent Coverage', 'E-value'  
    Naming convention ex: S204vsall.csv for the S204 genome  

2. A gene name .csv file which specifies information about each genome  
    Genename columns are 'Gene Name', 'locus_tag', 'Amino acid string'  
    Naming convention ex: S204_genes.csv for the S204 genome   
    The gene locus_tag should be the same format as the names appearing in the 'query' and 'subject' columns of the rawdata file 

## Output
#### Basic output
1. `core_genome.csv` - Contains rows of genes that shared across all submitted genomes
2. `core_names.csv` - Rows matching `core_genome.csv` that contain the gene names of the core genome
3. `all_aminos.csv` - A concatenated amino acid string of each shared gene from each genome
4. `{genome}_aminos.csv` - Rows matching `core_genome.csv` with the unique amino acid strings from each genome

#### Optional
1. `values_list.csv` - Output of all calculated values from geneparser runs, can be plotted
2. `core_pan.csv` - Output of core and pan genome sizes when incrementally adding genomes into the set, can be plotted
3. `pan_genome.csv` - Matrix of presence for the pan genome, per gene
4. `pan_genome_concat.csv` - Concat of pan genome presence for each organism

## Usage
The most basic run of geneparser.py is simple, just specify the input directory with -d  
`python3.6 geneparser.py -d /location/of/raw/data`  
This will calculate the shared genome at the default cutoffs (90/90/.0001) and place the results in `/location/of/raw/data/output`

#### Flags
-h, --help                  - show this help message and exit
-i, --percent-identity      - Percent identity cutoff for core genome
-c, --percent-coverage      - Percent coverage cutoff for core genome
-e, --expected-value        - Expected value cutoff for core genome
-p, --pan-genome            - Construct pan genome as part of output
-v, --verbose               - Increase output verbosity
-s, --save-gp               - Save gp file for test runs (takes a while)
-g, --core-pan-progression  - Generate csv of core and pan genome progression sizes
-u, --unique                - Generate list of uniques
-l, --list-values           - Append all found values w/ cutoffs to csv for graphing
-w, --overwrite-list        - Overwrite a previous values list
-x, --core-pan-intersection - Calculate intersection of core and pan genome
-r, --rule-out-similar      - Don't rule out similar genes from the same genome
-o, --output                - Desired output folder location
-d, --input-directory       - Folder containing raw data, output will be placed here as well if -o is not passed

## Utils

Extra scripts to expand usage and functionality of Geneparser

### genbank_build_blast.py

Creates the raw data files taken in by geneparser.py from raw `.gb` and `.gbk` files

### gp_plot.py

Creates some nice plots/graphs with the output of geneparser

 
