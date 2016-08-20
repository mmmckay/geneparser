# geneparser
Gene parsing tool for finding shared genes across genomes  
Requires python3, numpy, pandas  

Takes 2 files per genome, which can be generated using build_blast.py (experimental):  

1) A rawdata .csv file containing a single genome BLASTed against all other desired genome  
Rawdata columns are 'query', 'subject', 'Percent Identity', 'Percent Coverage', 'E-value'  
Naming convention ex: S204vsAll.csv for the S204 genome  

2) A gene name .csv file which specifies information about each genome  
Genename columns are 'Gene Name', 'Gene db_xref', 'Amino acid string for gene'  
Naming convention ex: S204_genes.csv for the S204 genome   

The gene db_xref should be the same format as the names appearing in the 'query' and 'subject' columns of the rawdata file  

Run with:  
$ python3 geneparser.py /User/raw 75 100 .0001  

Run in directory of raw folder to sample use  

Will run in current directory if no directory is specified, the numbers specify percent identity, percent coverage and e-value cutoffs   

Additonal commands  
-unw Specifies that the user has included only some of the files from the vsAll BLAST  
-v Verbose command, incomplete  
-pof Creates an output file that appends PI, PC and shared genes after each run, useful if you need to look at shared gene variation over multiple runs  
-t select the number of threads for the analysis to run with (experimental)  
