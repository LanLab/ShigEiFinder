# ShigEiFinder manuscript scripts

###clade_specific_gene_comninations.py
####This script takes in 4 files/folders:
>````commandline
>--presence roary_output.csv
>````
>the presence absence output csv from roary 

>````commandline
>--assign tab_delimited_file.txt
>````
>tab delimited file with 2 columns (genome name and cluster assignment) and a header line.

>````commandline
>--blastdb blastdb_path
>````
>Blast database including all genomes used in roary input

>````commandline
>--genomefolder folder_path
>````
>Folder containing one fasta file for each genome containing all genes for that genome (fasta headers must match roary gene names)

#### other settings:
````commandline
--falseneg FALSENEG   percentage false negatives allowed for any individual gene(rounds up) (default: 0)
-c CLADE, --clade CLADE
                    clade for which the gene/s are being selected (default: None)
--maxfp MAXFP         maximum percentage of false positive to allow for any individual gene (default: 10)
--combsize COMBSIZE   maximum number of loci to use in combination to get specific sets (default: 2)
--reportno REPORTNO   number of loci sets to report (default: 10)
--minalign MINALIGN   minimum alignment length to report BLAST result (default: 50)
--maxcombfp MAXCOMBFP
                    maximum total number of false positives for a set of genes (default: 20)
--cpus CPUS           number of threads to use (for BLAST) (default: 2)
--mingenesize MINGENESIZE
                    minimum length for a gene to be used (bp) (default: 200)
-o OUTPUT, --output OUTPUT
                    output file prefix (default: None)
````

####dependencies:
 - pandas
 - biopython


###prokka_genome_gene_from_roary.py

####This script takes in 4 files/folders:
>````commandline
>--presence roary_output.csv
>````
>the presence absence output csv from roary 

>````commandline
>--clusterassigns tab_delimited_file.txt
>````
>tab delimited file with 2 columns (genome name and cluster assignment) and a header line.

>````commandline
>--specificgenes tab_delimited_file.txt
>````
>tab delimited file with 2 columns (cluster name and specific gene roary orthologue group ID) and a header line.

>````commandline
>--genomefolder folder_path
>````
>Folder containing one fasta file for each genome containing all genes for that genome (fasta headers must match roary gene names)

#### other settings:
````commandline
--outprefix OUTPUT_PREFIX   prefix for output files
````

####dependencies:
 - biopython