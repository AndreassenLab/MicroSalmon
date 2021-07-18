# MicroSalmon

The MicroSalmon repository contains compresensive miRNA target and cis-regulatory motif prediction results for the Atlantic Salmon transcriptome, along with search scripts for exploring the dataset. Please cite Ramberg et al (DOI: To be added) if you make use of these data or programs in academic research.

## Dependencies

The included search scripts have no external dependencies beyond Python 3. Additionally, all the scripts are set up to be able to take arguments interactively simply by running the programs with no arguments, allowing for use without a command line terminal. All scripts must be in the same folder as the DATA folder to function.

## mRNA_Search.py
### Usage

```python mRNA_Search.py [-h] [-q QUERY] [-i QUERY_FILE] [-p PREFIX] [-s] [-c COMPLEXITY] [-u] [-a] [-m] [-t] [-r] ```  
Alternatively, arguments can be given interactively by running in script and writing the arguments in the provided field as you would for command line.
All outputs are stored in the OUTPUT folder.

### Arguments:
-h, --help Shows help message and exits

Search terms can be provided as an argument, or by giving a file with one search term per line, or both, in which case all terms from both inputs will be used:  
-q QUERY, --query QUERY                Query IDs, one or more accession numbers or SeqIDs separated by semicolons, e.g -q GIYK01000001;GIYK01000002 or -q SS1.1;SS1.2.  
-i QUERY_FILE, --query_file QUERY_FILE Query filename. Must be a list of accession numbers or SeqIDs, one per line.

Misc parameters:  
-p PREFIX, --prefix                    PREFIX Optional output prefix. If a prefix is provided, outputs will be stored using the name prefix.query.txt, otherwise they will simply be stored as query.txt.  
-s, --seqid                            Must be provided if the search term is a SeqID. Otherwise, the script will assume inputs are Accession numbers.  
-c COMPLEXITY, --complexity COMPLEXITY A demical number indicating minimum Trifonov Linguistic Complexity for Teiresias motifs to be included in the output. If not provided, the program defaults to 0.27.  

Output filters. The following arguments can be provided to selectively remove sections of the output file, as described in the section below.  
-u, --utr              Exclude complete 3'UTR sequence from output.  
-a, --annotation       Exclude gene and GO annotation from output.  
-m, --utrscan_motifs   Exclude UTRscan results from output.  
-t, --teiresias_motifs Exclude Teiresias results from output.  
-r, --mirna            Exclude miRNA target prediction results from output.
  
### Output Structure:

Header, lists the Accession number and SeqID of the input transcript, as well as the length of the 3'UTR.
```
Accession Number: GIYK01000002
SeqID in TSA: SS9.1
3'UTR length: 226
```

Complete 3'UTR sequence, can be removed from output using the argument -u.
```
3'UTR Sequence: CAACCAGGGTTAA...
```

Annotation, lists the gene and transcript annotation given by SQANTI, followed by the the gene symbol and description generated from the Blast analysis performed by Omicsbox, followed by the GO terms and IDs from the OmicsBox annotation. GO terms and IDs are provived as comma-separated lists, ordered with corresponding terms and IDs having the same position in their respective lists. Annotation can be removed from the output using the argument -a.
```
SQANTI Annotation Gene Symbol: LOC100316613
SQANTI Annotation Transcript ID: XM_014131490.1

BLAST Annotation Gene Symbol: mavs
BLAST Annotation Description: interferon promoter stimulating protein 1 isoform X1

GO Terms: activation of innate immune response,mitochondrion,...
GO IDs: GO:0002218,GO:0005739,...
```

UTRScan predictions, a header, line followed by all motifs predicted by UTRscan on separate lines. Each line contains the class of the motif, followed by the start and end positions of the motif in the 3'UTR, followed by the sequence of the motif. UTRscan predictions can be removed from the output using the argument -m.
```
UTRScan Motifs [Position in 3'UTR] Sequence:
K-BOX [1573,1580] CTGTGATG
K-BOX [3563,3570] CTGTGATC
```

Teiresias predictions, a header line, followed by all overrepresented motifs predicted by Teiresias passing the complexity filter. Each line contains the sequence of the motif, followed by the start and end positions of the motif in the 3'UTR, followed by the calculated linguistic complexity of the motif. Teiresias predictions can be removed from the output using the argument -t.
```
Teiresias Motifs [Position in 3'UTR] Linguistic Complexity(Threshold 0.27):
AAATGGCAC [1639,1648] 0.875
AATGGCACC [1640,1649] 1.0
```

miRNA target prediction, consists of a line listing the total number of miRNAs predicted to target the query transcript, followed by a list of all the miRNAs. For each predicted miRNA response element (MRE), the output then lists the number of the entry, followed by the name of the targeting miRNA and its mature sequence, folloqws by a list of prediction tools supporting the MRE prediction. After this, the output lists the minimum free energy for the miRNA/mRNA binding predicted by RNAhybrid, the position of the MRE in the 3'UTR, and a four line ascii illustration if the predicted miRNA/mRNA binding. The miRNA target prediction can be removed from the output by using the argument -r
```
Number of predicted MREs: 5
miRNAs predicted to bind MREs: ssa-miR-16bc-3p;ssa-miR-22a-3p;ssa-miR-22b-3p;ssa-miR-26a-4-3p;ssa-miR-26b-3p

MRE 1

miRNA: ssa-miR-16bc-3p
Mature miRNA sequence: CCCAATATTAGCTGTGCTGCTTC
Target prediction supported by: PITA, miRanda, TargetSpy, RNAhybrid

RNAhybrid output:
mfe: -24.9 kcal/mol

Position in 3'UTR: 180
Target 5' A                     A 3'
           AAGUAGUGCA CU AUAUUGG 
           UUCGUCGUGU GA UAUAACC 
miRNA  3' C          C  U       C 5'

MRE 2

miRNA: ssa-miR-22a-3p...
```

## miRNA_Search.py
### Usage
