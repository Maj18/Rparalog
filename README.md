# Rparalog:  
## Overview: 
Paralogs are homologous DNA/protein sequences that are found within species, they are the outcomes of gene duplicateion events. Gene duplications are important sources of genetic novelty and paralogs can also cause diseases (such as cancer). Rparalog is one program that identifies paralogs from an assembled genome based on protein sequences similarity using a relaxed reciprocal-best-hit BLAST strategy.

### Features: 
* A stringent **E-value** threshold (user can provided their own choice via a command line argument _-e_) plus the blast reciprocal-hit rule for accepting homologs (See Analysis 3) should have a good control of the prediction confidence. 
* Through command line argument (_-e_, see Usage), the user can decide the E-value threshold based on their study system.
* The paralog **annotation** results can also help the users to inspect their results.
* The protein sequences for all the members within a cluster will be collected into single files, to facilitate any follow-up analyses.
 
### Usage
* The program should be run as 
```bash
./src/Rparalog.py -g genomefile -f gtffile -p blastpfile -b namebase -e evalue & disown`
```
* Here, 
	* genomefile are an assemble genome sequence fasta file,
	* gtffile is a .gtf that is generated from a GeneMark gene prediction analysis, 
	* blastpfile is a .blastp file that comes from a blastp of the predicted protein sequence files again SwissProt database (uniprot database is also fine).
	* namebase for defining the output file names, 
	* evalue = the E-value of BLAST, evalue should be smaller than 1e-10 (this the e-value that is used for the self blast within the program). 

* e.g. 
```bash
./src/Rparalog.py -g data/Plasmodium_knowlesi.genome -f data/Plasmodium_knowlesi.gtf -p ./blastp/pk.blastp -b pk -e 1e-50 & disown
```
* (NBS. Please use the example command line for testing as it is given, because some of the steps (like the gffparse, makeblastdb, and blastp steps) in the program take too long time, therefore for testing, the output for those time-consuming steps have been included in the repository and therefore can be used directly for the following steps.)

## Software version:
* Python 3.7.3 (with sqlite3, argparse, subprocess installed)
* gffParse.pl: stands in bin folder, got from teacher
* makeblastdb 2.7.1+
* blastp 2.7.1+
* sqlite 3.27.2

## Analysis
1. Self BLASTP:
	A. Parse the .gtf file with the help of genome file using gffParse.p, a .faa (protein sequence fasta) file will be generated.
	B. Build a database from the .faa file.
	C. Blast the .faa file against its own database (from 1B).
2. Parse  the self .blastp file:
	* In this step, the .blastp file from 1C will be parsed to withdraw all the hits for each query, The output file cotains query-hit pairs.
3. Get reciprocal hits:
	* Only protein pairs that are each other hits within the 1C blastp analysis are kept as paralogous pairs.
4. Paralog clustering:
	* Paralogous pairs were clustered together into a paralog cluster if they have any shared members.
5. Annotate the predicted paralogs:
	1. parse the provide .blastP file (got via command-line arguments, it is from blast of the 1A .faa file against a SwissProt database ) and get the uniprot accession no. of the best hits for each query.
	2. Use the uniprot id from 5B to withdraw the corresponding pfam domains information and protein descriptions from the self-built databse SwissProti.sqlite (provided).
	3. Withdraw protein sequences from the 1A .faa file for all the members within a paralog cluster into one file and put in the paralog_seq folder.
	4. The SwissProt.sqlite database is built as follows:
		1. Download SwissProt database: `wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/ uniprot_sprot.dat.gz`
		2. Put each SwissProt record on one line (sprot.tab), createUniProtTab.pl -i uniprot_sprot.dat -o sprot
		3. Create the database: 
		```bash
		sqlite3 SwissProt.sqlite
		```
			.separator "\t"
			CREATE TABLE swissprot (accnr CHAR(7) PRIMARY KEY, description VARCHAR(240), taxid INTEGER(8), location VARCHAR(50), interpro VARCHAR(310), pfam VARCHAR(130), go_c VARCHAR(240), go_f VARCHAR(200), go_p VARCHAR(1100), ec VARCHAR(170), sequence TEXT);`
			.import sprot.tab swissprot`
			CREATE INDEX accnr ON swissprot (accnr);
			.quit
			
## Output files:
Among all the generated folders and files, three deserve special attention:
1. pk_paralog.out

	The_number_of_cluster_members |  Paralogous_copy1  |  Paralogous_copy2  |  .......
	---  |  ---  |  ---
	
2. pk_pannotate.out.

Gene_id  |  UniprotID_of_blastp_best_hit  |  PfamID_of_blastp_best_hit |  Functional_description_of_blastp_best_hit
---  |  ---  |  ---  |  ---  |

![alt text](https://github.com/Maj18/Rparalog/blob/master/pk_pannotate.out.png) 

3. paralog_seq folder: here the protein sequences of all the paralogous copies for each cluster will be collected into one single file.
	
			
		
	