# ShigEiFinder

This is a tool that is used to identify differentiate *Shigella*/EIEC 
using cluster-specific genes and identify the serotype using O-antigen/H-antigen genes. 
This pipeline can serotype over 59 *Shigella* and 22 EIEC serotypes using either assembled whole genomes 
or Whole Genome Sequencing (WGS) reads. The results are output in a tabular format which if saved as a 
file can be opened in Excel or other tabular programs. Also available as an [online tool](https://mgtdb.unsw.edu.au/ShigEiFinder/) and published in [microbial genomics](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000704).

Example:
````
#SAMPLE	ipaH	VIRULENCE_PLASMID	CLUSTER	SEROTYPE	O_ANTIGEN	H_ANTIGEN	NOTES
ERR1000679	+	1	CSD1	SD1	SD1		
````

---
# Installation 
## Dependencies
1. samtools (v1.10)
2. Python (v3.7.3 or above)
3. bwa (v0.7.17-r1188)
4. BLAST+ (v2.9.0)

## Option 1: Clone repository from gitHub
````
git clone https://github.com/LanLab/ShigEiFinder.git

cd ShigEiFinder

python setup.py install
````
Make sure that you have the dependencies installed.

## Option 2: Conda installation  
````
conda install -c conda-forge -c bioconda shigeifinder
````
[![Anaconda-Server Badge](https://anaconda.org/bioconda/shigeifinder/badges/installer/conda.svg)](https://conda.anaconda.org/bioconda) [![Anaconda-Server Badge](https://anaconda.org/bioconda/shigeifinder/badges/downloads.svg)](https://anaconda.org/bioconda/shigeifinder) [![Anaconda-Server Badge](https://anaconda.org/bioconda/shigeifinder/badges/version.svg)](https://anaconda.org/bioconda/shigeifinder)
---
# Usage
For genomes:
````
shigeifinder -i <inputs>...
````
For raw reads:
````
shigeifinder -r -i <read1> <read2>...
````

## Parameters Description
- ````-i```` : Input files list provided with their paths. Either genomes or read files can be used
- ````-r```` : Used to indicate that raw read files are input. Make sure that the reads are put in the order of read1 and then read2
- ````--hits```` : provides the genes set that was used to identify the cluster and serotype as well as the original BLAST/mapping results.
- ````--dratio````: displays the depth ratio of the depth of the cluster genes to the average depth of 7 HK genes
- ````--update_db````: updating the intermediate files for the genes database when new gene sequences have been added to the FASTA file
- ````-t````: number of threads used. The default is 4.
- ````--output````: output file to write to (if not used writes to stdout)
- ````--check````: check that dependencies are found in path
- ````--o_depth````: When using reads as input the minimum depth percentage relative to genome average for positive O antigen gene call (default 1)
- ````--ipaH_depth````: When using reads as input the minimum depth percentage relative to genome average for positive ipaH gene call (default 1)
- ````--depth````: When using reads as input the minimum read depth for non ipaH/Oantigen gene to be called (default 10)

---
# Troubleshooting
- If there is an error where there are no genes being mapped, there might be not enough memory given to the machine. Try giving more memory.
