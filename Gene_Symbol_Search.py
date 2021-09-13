#!/usr/bin/env python3

import os,sys
from argparse import ArgumentParser



parser = ArgumentParser()
parser.add_argument("-q", "--query", help="Query IDs, one Gene symbol, or a list of multiple Gene symbols in quotes separated by spaces, e.g -q \"ARGI2 mavs\"")
parser.add_argument("-i", "--query_file", help="Query filename. Must be a list of Gene Symbols, one per line. At least one of Input_File or Query is required.")
parser.add_argument("-p", "--prefix", help="Optional output prefix. If not included, output filename will be based only on the query ID.")
parser.add_argument("-c", "--complexity", default=0.27, type=float, help="Minimum Trifonov Linguistic Complexity for Teiresias motifs. (Default: 0.27)")
parser.add_argument("-u", "--utr", default=False, action="store_true", help="Exclude complete 3'UTR sequence from output.")
parser.add_argument("-a", "--annotation", default=False, action="store_true", help="Exclude gene and GO annotation from output.")
parser.add_argument("-m", "--utrscan_motifs", default=False, action="store_true", help="Exclude UTRscan results from output.")
parser.add_argument("-t", "--teiresias_motifs", default=False, action="store_true", help="Exclude Teiresias results from output.")
parser.add_argument("-s", "--target_summary", default=False, action="store_true", help="Exclude summary on target transcripts, genes and GOs from output.")
parser.add_argument("-r", "--mirna", default=False, action="store_true", help="Exclude miRNA target prediction results from output.")

if not len(sys.argv) > 1:
    print("No command line arguments provided. Please provide arguments below. For help, type -h")
    sys.argv.extend(input("Args: ").split())
    
args = parser.parse_args()

if not args.query and not args.query_file:
    print("No query provided. Include either a query or a query file.")
    sys.exit(-1)

if args.query_file:
    if not os.path.exists(args.query_file):
        print("Query file does not exist")
        sys.exit(-1)
f=open("DATA/SQANTI_OmicsBox_Annotation.tsv","r")
lines = f.readlines()
f.close()

annotation_dict = {}
seqid_to_accs = {}
gene_dict = {}
for line in lines[1:]:
    annotation_dict[line.split("\t")[1]]=line.split("\t")
    seqid_to_accs[line.split("\t")[1]]=line.split("\t")[0]
    if "-NA-" not in line.split("\t")[4] and line.split("\t")[4] != "":
        if line.split("\t")[4] not in gene_dict:
            gene_dict[line.split("\t")[4]] = set([line.split("\t")[1]])
        else:
            gene_dict[line.split("\t")[4]].add(line.split("\t")[1])
    if "novelGene" not in line.split("\t")[2] and line.split("\t")[2] != "" and "CG" not in line.split("\t")[1]:
        if line.split("\t")[2] not in gene_dict:
            gene_dict[line.split("\t")[2]] = set([line.split("\t")[1]])
        else:
            gene_dict[line.split("\t")[2]].add(line.split("\t")[1])

f=open("DATA/uscan_output.txt","r")
lines = f.readlines()
f.close()

uscan_dict = {}

for line in lines[48:]:
    if line.split(" ")[0].split(";")[0] not in uscan_dict:
        for ID in line.split(" ")[0].split(";"):
            uscan_dict[ID] = [line.split(" : ")[1] + " " + line.split(" : ")[2]]
    else:
        for ID in line.split(" ")[0].split(";"):
            uscan_dict[ID] += [line.split(" : ")[1] + " " + line.split(" : ")[2]]

f=open("DATA/mRNA_3UTR.fasta","r")
lines = f.readlines()
f.close()

utr_dict = {}

for line in lines:
    if line[0] == ">":
        ID = line.split(";")[1].strip()
    else:
        utr_dict[ID] = line.strip()

f=open("DATA/CD-Hit_Clusters.txt","r")
lines = f.readlines()
f.close()

counter = 0
tnum_dict = {}

for line in lines:
    for ID in line.strip("\n").split(";"):
        tnum_dict[ID] = str(counter)
    counter += 1

f=open("DATA/Teiresias_k1000_prob_cutoff_5.txt","r")
lines = f.readlines()
f.close()

teiresias_dict = {}

for line in lines:
    s = line.split("\t")[2].split(" ")[0]
    c = 0
    for base in ["A","T","G","C"]:
        if base in s:
            c += 1
    c = c/4
    for k in range(2,len(s)):
        words = set()
        for pos in range(0,len(s)-k+1):
            words.add(s[pos:pos+k])
        c = c*len(words)/min(4**k,len(s)-k+1)
    if c > args.complexity:
        counter = 0
        for numb in line.strip("\n").split("\t")[2].split(" ")[1:]:
            if counter == 0:
                ID = numb
                counter = 1
            else:
                if ID not in teiresias_dict:
                    teiresias_dict[ID] = [s + " [" + str(int(numb)+1) + "," + str(int(numb)+1+len(s)) + "] " + str(c) + "\n"]
                else:
                    teiresias_dict[ID] += [s + " [" + str(int(numb)+1) + "," + str(int(numb)+1+len(s)) + "] " + str(c) + "\n"]
                counter = 0


RNAhybrid_dict = {}
for numb in ["1","2","3","4"]:
    f=open("DATA/RNAhybrid_target_prediction_part_"+numb+".txt","r")
    lines = f.readlines()
    f.close()

    for line in lines:
        if line.split(":")[0].split(";")[0] not in RNAhybrid_dict:
            for ID in line.split(":")[0].split(";"):
                RNAhybrid_dict[ID] = [line.strip("\n").split(":")]
        else:
            for ID in line.split(":")[0].split(";"):
                RNAhybrid_dict[ID] += [line.strip("\n").split(":")]
                
f=open("DATA/RNAhybrid_plus_2.txt","r")
lines = f.readlines()
f.close()

prediction_dict = {}

for line in lines:
    for ID in line.split("\t")[1].split(";"):
        prediction_dict[ID+line.split("\t")[0]] = line.split("\t")[3]

f=open("DATA/miRNAome.fa","r")
lines = f.readlines()
f.close()

miRNA_dict = {}

for line in lines:
    if line[0] == ">":
        ID = line[1:].strip()
    else:
        miRNA_dict[ID] = line

queries = []

if args.query:
    queries += args.query.split(" ")

if args.query_file:
    f=open(args.query_file,"r")
    lines = f.readlines()
    f.close()
    for line in lines:
        queries.append(line.strip())

while "" in queries:
    queries.remove("")
        
for counter, query in enumerate(queries):
    print("Searching query " + str(counter+1) + " of " + str(len(queries)))
    skip = 0
    if query in gene_dict:
        print(query)
    else:
        print(query + " not found. Check spelling.")
        skip = 1
    if skip == 0:
        output = "Input Gene Symbol: " + query + "\n"
        if not args.target_summary:
            accs = []
            GOIDs = []
            GOTs = []
            MREs = set()
            for seqid in sorted(list(gene_dict[query])):
                accs.append(seqid_to_accs[seqid])
                if seqid in RNAhybrid_dict:
                    for MRE in RNAhybrid_dict[seqid]:
                        MREs.add(MRE[2])
                for GOID in annotation_dict[seqid][6].split(","):
                    if GOID not in GOIDs:
                        GOIDs.append(GOID)
                for GOT in annotation_dict[seqid][7].split(","):
                    if GOT not in GOTs:
                        GOTs.append(GOT)
            MREs = sorted(list(MREs))
            output += "\nNumber of Annotated Transcripts: " + str(len(gene_dict[query])) + "\nTranscript Accessions: " + ";".join(accs) + "\nTranscript SeqIDs in TSA: " + ";".join(sorted(list(gene_dict[query])))
            output += "\n\nNumber of predicted MREs for Gene: " + str(len(MREs)) + "\nmiRNAs predicted to bind MREs: " + ";".join(MREs) + "\n\nNumber of Unique Annotated GO Terms: " + str(len(GOTs)) + "\nGO Terms: " + ";".join(GOTs) + "\nGO IDs: " + ";".join(GOIDs)+"\n"
        for counter, seqid in enumerate(sorted(list(gene_dict[query]))):
            output += "\nTranscript " + str(counter+1) + ":\n\nAccession Number; " + seqid_to_accs[seqid] + "\nSeqID in TSA: " + seqid + "\n3'UTR length: " + str(len(utr_dict[seqid])) +  "\n"
            if not args.utr:
                output += "\n3'UTR Sequence: " + utr_dict[seqid] + "\n"
            if not args.annotation:
                output += "\nSQANTI Annotation Gene Symbol: " + annotation_dict[seqid][2] + "\nSQANTI Annotation Transcript ID: " + annotation_dict[seqid][3] + "\n\nBLAST Annotation Gene Symbol: " + annotation_dict[seqid][4] + "\nBLAST Annotation Description: " + annotation_dict[seqid][5] + "\n\nGO Terms: " + annotation_dict[seqid][7] + "\nGO IDs: " + annotation_dict[seqid][6] +"\n"
            if not args.utrscan_motifs:
                output += "\nUTRScan Motifs [Position in 3'UTR] Sequence:\n"
                if seqid in uscan_dict:
                    motifs = []
                    for motif in uscan_dict[seqid]:
                        motifs.append(motif)
                    motifs.sort()
                    for motif in motifs:
                        output += motif
                else:
                    output += "No hits\n"
            if not args.teiresias_motifs:
                output += "\nTeiresias Motifs [Position in 3'UTR] Linguistic Complexity(Threshold " + str(args.complexity) + "):\n"
                if tnum_dict[seqid] in teiresias_dict:
                    motifs = []
                    for motif in teiresias_dict[tnum_dict[seqid]]:
                        motifs.append(motif)
                    motifs.sort()
                    for motif in motifs:
                        output += motif
                else:
                    output += "No hits\n"
            if not args.mirna:
                if seqid in RNAhybrid_dict:
                    output += "\nNumber of predicted MREs for " + seqid + ": " + str(len(RNAhybrid_dict[seqid])) + "\nmiRNAs predicted to bind MREs: "
                    mirnas = set()
                    for mirna in RNAhybrid_dict[seqid]:
                        mirnas.add(mirna[2])
                    mirnas = list(mirnas)
                    mirnas.sort()
                    output += ";".join(mirnas) + "\n"
                    for counter, mirna in enumerate(RNAhybrid_dict[seqid]):
                        output += "\nMRE " + str(counter+1) + "\n\nmiRNA: " + mirna[2] + "\nMature miRNA sequence: " + miRNA_dict[mirna[2]] + "Target prediction supported by: " + prediction_dict[seqid+mirna[2]]
                        output += "\nRNAhybrid output:\nmfe: " + mirna[4] + " kcal/mol\n\nPosition in 3'UTR: "  + mirna[6] + "\nTarget 5' " + mirna[7]
                        output += " 3'\n          " + mirna[8] + "\n          " + mirna[9] + "\nmiRNA  3' " + mirna[10] + " 5'\n"
                    
            else:
                output += "\nNo MREs\n"



        if args.prefix:
            outputfile = "OUTPUT/" + args.prefix + "." + query + ".txt"
        else:
            outputfile = "OUTPUT/" + query + ".txt"
        f=open(outputfile,"w")
        f.write(output)
        f.close()
if args.prefix:
    print("Search complete. Results are found in the OUTPUT folder with the prefix " + args.prefix + ".")
else:
    print("Search complete. Results are found in the OUTPUT folder.")
