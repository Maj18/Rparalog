# Rparalog Manual
This manual corresponds to version 1.1

## Introduction
Paralogs are homologous DNA/protein sequences that are found within species, they are the outcomes of gene duplicateion events. Gene duplications are important sources of genetic novelty and paralogs can also cause diseases (such as cancer). Rparalog is one program that identifies paralogs from an assembled genome based on protein sequences similarity using a relaxed reciprocal-best-hit BLAST strategy.

### Features: 
* A stringent **E-value** threshold plus the blast reciprocal-hit rule for accepting homologs (See Analysis 3) should have a good control of the prediction confidence. 
* Through command line argument (_-e_, see Usage), the user can decide the E-value threshold based on their study system.
* The paralog **annotation** results can also help the users to inspect their results.
* The protein sequences for all the members within a cluster will be collected into single files, to facilitate any follow-up analyses.

## Installation
### Prerequisites
To run Rparalog, you need:
* Python 3.7.3 (with sqlite3, argparse, subprocess installed)
* makeblastdb 2.7.1+
* blastp 2.7.1+
* sqlite 3.27.2
* a simple SwissProt.sqlite database should be built beforehand (see below for more instruction), and SwissProt.sqlite should be put in the folder src. 
		1. Download SwissProt database: `wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/ uniprot_sprot.dat.gz`
		2. Put each SwissProt record on one line (sprot.tab), 
		```bash
		createUniProtTab.pl -i uniprot_sprot.dat -o sprot
		```
		3. Create the database: 
		```bash
		sqlite3 SwissProt.sqlite
		```
			.separator "\t"
			CREATE TABLE swissprot (accnr CHAR(7) PRIMARY KEY, description VARCHAR(240), taxid INTEGER(8), location VARCHAR(50), interpro VARCHAR(310), pfam VARCHAR(130), go_c VARCHAR(240), go_f VARCHAR(200), go_p VARCHAR(1100), ec VARCHAR(170), sequence TEXT);
			.import sprot.tab swissprot
			UPDATE swissprot SET location = null WHERE location = "null";
			UPDATE swissprot SET interpro = null WHERE interpro = "null";
			UPDATE swissprot SET pfam = null WHERE pfam = "null";
			UPDATE swissprot SET go_c = null WHERE go_c = "null";
			UPDATE swissprot SET go_f = null WHERE go_f = "null";
			UPDATE swissprot SET go_p = null WHERE go_p = "null";
			UPDATE swissprot SET ec = null WHERE ec = "null";
			CREATE INDEX accnr ON swissprot (accnr);
			.quit
   * (In databases a missing value is represented by the null value. However when we import data to a SQLite table it’s not possible to set a value to be absent. In our data we have represented missing value as the string null. We have to convert these values to true null values, that's why we did the UPDATE steps above)

### Building and installing 
* If you have git installed on your computer, you can simply get the software by 
```bash
git clone https://github.com/Maj18/Rparalog.git
```
* Otherwise, you can download the software folder from https://github.com/Maj18/Rparalog.
* Change directory into the extracted folder e.g. via cd Rparalog or cd Rparalog-master (if downloaded directly)
* you can now run ./src/Rparalog.py ... directly.

## Quick Start
./src/Rparalog.py -p proteinfile  -b namebase -e evalue & disown

## Quick Tutorial
Rparalog assumes that you have all your protein sequences in FASTA format, the software package contains an example, namely pk.faa. The full command line for the example files would thus be `./src/Rparalog.py -p pk.faa -b pk -e 1e-50 & disown`. (Please be aware that examples of SwissProt.sqlite database, .blastp (against SwissProt database) and self-database and .blast (against self-built database) have been provided for quick test, due to the uploading size limitation of Github, some of these files are not complete, therefore the ..._pannotate.out will have lots of missing data, just you know it.

Below is a step by step explanation of how the software works:

### Analysis
**Part I.**
1. Blast the user-provide protein sequence file (.faa) against Swissprot database, the outcome will be used for the paralog annotation later on.
2. Self BLASTP:
	A. Build a database from the .faa file.
	B. Blast the .faa file against its own database (from 2A).

**PartII**
3. Parse  the self .blastp file:
	* In this step, the .blastp file from 2B will be parsed to withdraw all the hits for each query, The output file contains query-hit pairs.
4. Get paralogous pairs:
	* Only protein pairs that are each other's hit within the 2B blastp analysis are kept as paralogous pairs.
5. Paralog clustering:
	* Paralogous pairs were clustered together into paralog clusters if they have any shared members.

**PartIII**
6. Annotate the predicted paralogs:
	A. parse the .blastP file (from step 1) and get the uniprot accession no. of the best hit for each query.
	B. Use the uniprot id from 6A to withdraw the corresponding pfam domain and protein description information from the self-built databse SwissProti.sqlite (see **Prerequisites**).
	C. Withdraw protein sequences from the  .faa file for all the members within a paralog cluster into one file and put it in the paralog_seq folder.
			
### Output files:
Among all the generated folders and files, three deserve special attention:
1. {namebase}_paralogCluster.out

	The_number_of_cluster_members |  Paralogous_copy1  |  Paralogous_copy2  
	---  |  ---  |  ---
	
2. {namebase}_pannotate.out.

Gene_id  |  UniprotID_of_blastp_best_hit  |  PfamID_of_blastp_best_hit |  Functional_description_of_blastp_best_hit
---  |  ---  |  ---  |  ---  |

3. paralog_seq folder: here the protein sequences of all the paralogous copies for each cluster will be collected into one single file.


## Usage
* The program should be run as 
```bash
./src/Rparalog.py -p proteinfile  -b namebase -e evalue & disown
```
* Here, 
-p | proteinpfile | a .faa file that include all the protein sequences of a genome.
-b | namebase |  prefix for all output file names, 
-e | evalue |  = the E-value of BLAST, evalue should be smaller than 1e-10. 
---  |  ---  |  ---  

## Caution!

It should be noted that, the Rparalog program can only predict what paralogs are there in a genome, but can not tell their exact copy number, the variation of which has been found to cause many different kinds of diseases. This is because nowadays, the assembled genomes are mostly haploid, due to the assembling challenge in separating between homologous chromosomes. However, as the fast development of long-fragment sequencing technology, it may be possible one day to assemble homologous chromosomes separately, then it would be much easier to use Rparalog to decide the copy numbers of paralogs with high confidence.			
		
	