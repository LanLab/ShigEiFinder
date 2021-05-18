# ShigEiFinder manuscript scripts

This folder contains scripts used in the shigeifinder [manuscript](https://www.biorxiv.org/content/10.1101/2021.01.30.428723v3) to identify and extract specific gene sets.

### clade_specific_gene_combinations.py
This script was used to identify specific gene sets for each cluster from the pan genomes of the identification dataset. The script ran on one cluster at a time. The script takes in 4 inputs, a roary presence absence file, a genome cluster assignment file, the genomes of all isolates, the annotated genes in all genomes (as used in roary). The script first identified individual candidate genes that were present in all isolates of the target cluster (true positives) and were present only in a percentage of non-target cluster isolates (false positives). For the list of candidates each combination of genes was tested to see whether all are found in the same false
positive strain. If a set of genes are never all found together then that set of genes is reported as a result. The size of the gene combinations starts at 1 for the whole list and increases progressively. At each size, successful sets of genes were reported until the total number of reported sets equals the maximum specified in the settings. Additionally, if a successful set of 2 genes (for example) was found within a subsequent set of 3 genes that three gene set was excluded because the additional gene provides no benefit.

This script takes in 4 files/folders:
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

other settings:
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

dependencies:
 - pandas
 - biopython


### prokka_genome_gene_from_roary.py

This script extracts specific gene set sequences for sets produced by clade_specific_gene_combinations.py. The script accepted 4 inputs: the presence absence roary output csv, the annotated genes in all genomes (as used in roary), a list of cluster specific genes sets and their corresponding cluster, a list of genome ids and their corresponding cluster. An output prefix is also required. The script will: 
 - select a representative genome from each cluster
 - identify the roary orthologue group that contains a given specific gene
 - retrieve the gene ID for that orthologue group and the representative genome
 - extract the gene from the genes fasta file for that genome
 - save the specific gene to an output file (output prefix)
 - produce a summary file of genes retrieved (output prefix)

This script takes in 4 files/folders:
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

other settings:
````commandline
--outprefix OUTPUT_PREFIX   prefix for output files
````

dependencies:
 - biopython
