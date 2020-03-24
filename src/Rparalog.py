#!/usr/bin/env python3

#In this version, the annotation of paralogs makes use of the genome annotation results (blastp of all the predicted protein sequences file against SwissProt database)

import os
import sys
import sqlite3
import argparse
from os import path
from subprocess import Popen, PIPE, STDOUT

if len(sys.argv) != 7:
    sys.exit("ERROR! The program should be run as ./src/Rparalog.py -p proteinfile -b namebase -e evalue & disown, \
        e.g. ./src/Rparalog.py -p pk.faa -b pk -e 1e-50 & disown")

#Parse the command line arguments using argparse
parser = argparse.ArgumentParser(description = 'This program will predict gene paralogs for a given set protein sequences using a relaxed reciprocal-best-hit blast strategy')
parser.add_argument(
    '-p',
    dest = "proteinfile",
    metavar = 'PROTEINFILE',
    type = str,
    default = True)
parser.add_argument(
    '-b',
    metavar = "NAMEBASE",
    dest = 'namebase',
    type = str,
    default = True)
parser.add_argument(
    '-e',
    metavar = "EVALUE",
    dest = 'evalue',
    type= float,
    default = True)
args = parser.parse_args()
proteinfile = args.proteinfile
namebase = args.namebase
evalue = args.evalue

def selfblastp(proteinfile,namebase):
    #BlastP the given proteins against SwissProt database, the output will be used at the paralog annotation step that will come later on.
    if not path.exists(f"blastp/{namebase}.blastp"):
        runCommand = f"blastp -query {proteinfile} -out blastp/{namebase}.blastp -db SwissProt"
        p1 = Popen(runCommand.split(), stdout=PIPE, close_fds=True, shell=False)
        p1.wait()    
    #make a databse from the given protein sequences
    if not path.exists(f"database/{namebase}.phr"):
        runCommand2 = f"makeblastdb -in {proteinfile} -out database/{namebase}Aa -dbtype prot"
        p2 = Popen(runCommand2.split(), stdout=PIPE, close_fds=True, shell=False)
        p2.wait()
    #Blastp the given protein sequences against its own datebase.
    if not path.exists(f"blastp/{namebase}_vs_{namebase}.blastp"):
        runCommand3 = f"blastp -query {proteinfile} -db database/{namebase}Aa -evalue 1e-10 -out blastp/{namebase}_vs_{namebase}.blastp -num_descriptions 50 -num_threads 2"
        p3 = Popen(runCommand3.split(), stdout=PIPE, close_fds=True, shell=False)
        p3.wait()
    

#Parse the self-blastp file acquired above to get all hits for each query sequence
#Please be aware that with the command line argument evalue, the user can have a control of the degree of sequence similarity between the query and target sequences
#Only hits with an e-value smaller than the provided evalue will be accepted for the the follow-up analyses.
def blastParser(namebase, evalue):
    with open (f"blastp/{namebase}_vs_{namebase}.blastp", "r") as fin, open (f"{namebase}_query_hit.out", "w") as fout:
        for line in fin:
            if "Query=" in line:
                hit = ""
                Query_line = line.rstrip().split()[1]
            elif line.startswith(">") and (Query_line != line.rstrip().split()[1]):
                hit = line.split()[1]
                EN = 0
            elif "Expect =" in line and hit:
                if EN==0:
                    EVAL = float(line.split(" ")[9][:-1])
                    if EVAL < evalue:
                        print(Query_line, hit, sep="\t", file=fout)
                        EN += 1
                        hit = ""


#Only protein pairs that are each other's blast hit (given the provided evalue above) will be accepted as paralogous pairs
def getParalogPair():
    query_hit = []
    with open(f"{namebase}_query_hit.out", "r") as fin2:
        for line2 in fin2:
            if line2:
                line2 = line2.rstrip().split("\t")
                query_hit.append(line2)
    duplicates = []
    with open (f"{namebase}_paralogPair.out", "w") as foutr:
        for [q, h] in query_hit:
            rev = [h, q]
            if rev in query_hit:
                if rev not in duplicates:
                    forw = [q,h]
                    forw.sort(key = lambda x: int(x.split("_")[0]))
                    duplicates.append(forw)
                    print(forw[0], forw[1], sep="\t", file=foutr)
    return duplicates


#We traverse the pairwise relationships (between the paralogous pairs) to find the maximally connected clusters that are disjoint from one another. 
#Paralogous pairs were clustered together if they have any shared members.
def paralogClustering():
    paralogs = []
    rest = paralogPair
    l = len(rest)
    while l>0:   
        paralog = set()
        paralog.add(rest[0][0])
        paralog.add(rest[0][1])
        del rest[0]  
        l = len(rest)
        if l == 0:
            break
        else:
            delete = []
            for i in range(l):
                if (rest[i][0] in paralog) or (rest[i][1] in paralog):
                    paralog.add(rest[i][0])
                    paralog.add(rest[i][1])
                    delete.append(i)
            delete.reverse()
            for n in delete:
                del rest[n]
            l = len(rest)
        paralogs.append(paralog)
    with open (f"{namebase}_paralogClusters.out","w") as fout2:
        for paralog in paralogs:
            print(len(paralog), file=fout2, end="\t")
            for copy in paralog:
                print(copy, end="\t", file=fout2)
            print(file=fout2)


#Here I provide the annotation information for the predicted paralogs for the users, it's up to them how to interpret the results
#Here the uniprot accession no., pfam id and description of the best blastp hit (against Swissprot database) will be collected for each predicted paralogs
#In order to do what I want, I have built my own SQLite database, called SwissProt.sqlite
#If time allows, it's better to the blastp against the Uniprot database, by doing that, one can get more gene annotated.
#First, parse the .blastP file (from blastp of the .faa file against Swissprot database) and get the uniprot accession no. of the best hits for each query
def blastpParserForAnnotate():
    with open (f"./blastp/{namebase}_besthit.swissprot", "w") as fout3, open(f"./blastp/{namebase}.blastp", "r") as blastpfile:
        for line in blastpfile:
            if "Query=" in line:
                Query_line = line.rstrip().split()[1]
                target=0
            elif (line.startswith(">")) and (target == 0):
                if Query_line != line.rstrip().split()[1]:
                    print(Query_line, file=fout3, end=" ")
                    print(line.split()[1].split("|")[1], file=fout3)
                    target += 1
                else: 
                    target=0

#The second annotation step is to use the uniprot id of the best hit from the blastp step to withdraw pfam and description information from the SQLite database
def pannotate(proteinfile):
    with open(f"./blastp/{namebase}_besthit.swissprot", "r") as fin3, open (f"{namebase}_paralogClusters.out", "r") as fin4, open(f"{namebase}_pannotate.out", "w") as fout4:
        #withdraw the uniprot id of the best hit from the blastp step
        uniprot = {}
        for line3 in fin3:
            line3 = line3.rstrip().split()
            uniprot[line3[0]] = line3[1]
        n = 1
        #withdraw pfam and description information from the SQLite database based on uniprot id
        for line4 in fin4: 
            line4 = line4.rstrip().split()
            print("#######Paralog_cluster", n, file=fout4)
            with open(f"paralog_seq/{namebase}_paralogCluster{n}.faa", "w+") as foutx: #print the protein sequences for each paralog cluster into a file (in the paralog_seq folder)
                for copy in line4[1:]:
                    runCommandx = f"grep -A 1 {copy} {proteinfile}"
                    Popen(runCommandx.split(), stdout=foutx, stderr=PIPE, universal_newlines=True).communicate()
                    if (copy not in uniprot) or (not uniprot[copy]):
                        print(copy, "null", sep="\t", file=fout4)
                    else:
                        print(copy, uniprot[copy], sep= "\t", end="\t", file=fout4)
                        uniprotID = (uniprot[copy],)
                        conn = sqlite3.connect('./src/SwissProt.sqlite') #Make a connection object
                        c = conn.cursor() #Create a cursor object
                        c.execute('SELECT pfam, description FROM swissprot WHERE accnr=?;', uniprotID)
                        row=""
                        for row in c:
                                print("\t".join(map(str, row)), end="\n", file=fout4)
                        if not row:
                            print(file=fout4)
            n += 1


if __name__ == "__main__":
    selfblastp(proteinfile, namebase)
    blastParser(namebase, evalue)
    paralogPair=getParalogPair()
    paralogClustering()
    blastpParserForAnnotate()
    pannotate(proteinfile)

