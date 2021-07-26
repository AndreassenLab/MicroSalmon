# MicroSalmon

The MicroSalmon repository contains compresensive miRNA target and cis-regulatory motif prediction results for the Atlantic Salmon transcriptome, along with search scripts for exploring the dataset. Please cite Ramberg et al (DOI: 10.3389/fgene.2021.656334 and to be added) if you make use of these data or programs in academic research, and include this citation notice with any modified version of this database or any other work including all or part of it.

## Table of contents
- [Dependencies](#dependencies)
- [General usage note](#general-usage-note)
- [mRNA_Search.py](#mrna_searchpy)
  * [Usage](#usage)
  * [Arguments:](#arguments)
  * [Output Structure:](#output-structure)
- [miRNA_Search.py](#mirna_searchpy)
  * [Usage](#usage-1)
  * [Arguments:](#arguments-1)
  * [Output Structure:](#output-structure-1)
- [Gene_Symbol_Search.py](#gene_symbol_searchpy)
  * [Usage](#usage-2)
  * [Arguments:](#arguments-2)
  * [Output Structure:](#output-structure-2)
- [GO_ID_Search.py](#go_id_searchpy)
  * [Usage](#usage-3)
  * [Arguments:](#arguments-3)
  * [Output Structure:](#output-structure-3)
- [DATA file structure](#data-file-structure)
  * [SQANTI_OmicsBox_Annotation.tsv](#sqanti_omicsbox_annotationtsv)
  * [miRNAome.fa](#mirnaomefa)
  * [mRNA_3UTR.fasta](#mrna_3utrfasta)
  * [CD-Hit_Clusters.txt](#cd-hit_clusterstxt)
  * [RNAhybrid_target_prediction_part_X](#rnahybrid_target_prediction_part_x)
  * [RNAhybrid_plus_2.txt](#rnahybrid_plus_2txt)
  * [Teiresias_k1000_prob_cutoff_5_no_miRNA_seed.txt](#teiresias_k1000_prob_cutoff_5_no_mirna_seedtxt)
  * [uscan_output.txt](#uscan_outputtxt)


## Dependencies

The included search scripts have no external dependencies beyond Python 3.  
Additionally, all the scripts are set up to be able to take arguments interactively simply by running the programs with no arguments, allowing for use without a command line terminal.  
All scripts must be in the same folder as the DATA and OUTPUT folders to function.

## General usage note
All miRNA, Accession Number, SeqID, gene symbol and GO ID queries are case sensitive, and must be written as shown in the files SQANTI_OmicsBox_Annotation.tsv or miRNAome.fa, or in any of the search outputs, for the search scripts to recognize them.

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

```python miRNA_Search.py [-h] [-q QUERY] [-i QUERY_FILE] [-p PREFIX] [-c COMPLEXITY] [-u] [-a] [-m] [-t] [-s] ```  
Alternatively, arguments can be given interactively by running in script and writing the arguments in the provided field as you would for command line.
All outputs are stored in the OUTPUT folder.

### Arguments:
-h, --help Shows help message and exits

Search terms can be provided as an argument, or by giving a file with one search term per line, or both, in which case all terms from both inputs will be used:  
-q QUERY, --query QUERY                Query IDs, one or more mature miRNA IDs separated by semicolons, e.g -q ssa-let-7a-4-5-3p;ssa-let-7a-3-5-5p  
-i QUERY_FILE, --query_file QUERY_FILE Query filename. Must be a list of mature miRNA IDs, one per line. At least one of Input_File or Query is required.  

Misc parameters:  
-p PREFIX, --prefix                    PREFIX Optional output prefix. If a prefix is provided, outputs will be stored using the name prefix.query.txt, otherwise they will simply be stored as query.txt.  
-c COMPLEXITY, --complexity COMPLEXITY A demical number indicating minimum Trifonov Linguistic Complexity for Teiresias motifs to be included in the output. If not provided, the program defaults to 0.27.  

Output filters. The following arguments can be provided to selectively remove sections of the output file, as described in the section below.  
-u, --utr              Exclude complete 3'UTR sequence from output.  
-a, --annotation       Exclude gene and GO annotation from output.  
-m, --utrscan_motifs   Exclude UTRscan results from output.  
-t, --teiresias_motifs Exclude Teiresias results from output.  
-s, --target_summary  Exclude summary on target transcripts, genes and GOs from output.  

### Output Structure:

Header, lists the name of the input miRNA, as well as the mature miRNA sequence.
```
Input miRNA: ssa-miR-10d-3p
Mature miRNA sequence: CAAATTCGCCTCTACGGGA
```

Summary data, the first subsection consists of two lines with the total number of target mRNA transcripts predicted for the miRNA, and the number of unique different target genes, followed by corresponding lists of Accession numbers, SeqIDs, and gene annotation for each transcript, in the same order. Gene annotation defaults to the OmicsBox annotation if present, followed by the SQANTI annotation, or lists --NA-- if no annotation was found. The second subsection shows the total number of unique GO terms associated with all the target miRNAs, and a non-redundant list of unique GO terms. Summary can be removed from the output using the flag -s.
```
Number of Target Transcripts: 196
Number of unique target genes: 82
Target Transcript Accessions: GIYK01066742;GIYK01067131;GIYK01067138;...
Target Transcript SeqIDs in TSA: CG132.2;CG806.1;CG821.2;...
Target Transcript Annotated Genes: --NA--;--NA--;LOC106588121;...

Number of unique GO Terms in target transcripts: 134
GO Terms: ubiquitin-protein transferase activity;aminoacyl-tRNA ligase activity;...
```
mRNA header, each individual mRNA entry begins with a header listing the number of the entry, followed by the Accession number, SeqID, and and length of the 3'UTR for the transcript in the entry.
```
Transcript 1:

Accession Number; GIYK01066742
SeqID in TSA: CG132.2
3'UTR length: 782
```
Complete 3'UTR sequence, can be removed from output using the argument -u.
```
3'UTR Sequence: ACCTACCCAC...
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

miRNA target prediction, consists of a list of prediction tools supporting the MRE prediction, followed by the minimum free energy for the miRNA/mRNA binding predicted by RNAhybrid, the position of the MRE in the 3'UTR, and a four line ascii illustration if the predicted miRNA/mRNA binding.
```
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

## Gene_Symbol_Search.py
### Usage

```python Gene_Symbol_Search.py [-h] [-q QUERY] [-i QUERY_FILE] [-p PREFIX] [-c COMPLEXITY] [-u] [-a] [-m] [-t] [-s] [-r] ```  
Alternatively, arguments can be given interactively by running in script and writing the arguments in the provided field as you would for command line.
All outputs are stored in the OUTPUT folder.

### Arguments:
-h, --help Shows help message and exits

Search terms can be provided as an argument, or by giving a file with one search term per line, or both, in which case all terms from both inputs will be used. gene symbols are case sensitive, and must be written as given in other search outputs or in the SQANTI_OmicsBox_Annotation.tsv file in the DATA folder:  
-q QUERY, --query QUERY                Query IDs, one or more Gene symbols separated by semicolons, e.g -q ARGI2;mavs  
-i QUERY_FILE, --query_file QUERY_FILE Query filename. Must be a list of Gene Symbols, one per line. At least one of Input_File or Query is required.  

Misc parameters:  
-p PREFIX, --prefix                    PREFIX Optional output prefix. If a prefix is provided, outputs will be stored using the name prefix.query.txt, otherwise they will simply be stored as query.txt.  
-c COMPLEXITY, --complexity COMPLEXITY A demical number indicating minimum Trifonov Linguistic Complexity for Teiresias motifs to be included in the output. If not provided, the program defaults to 0.27.  

Output filters. The following arguments can be provided to selectively remove sections of the output file, as described in the section below.  
-u, --utr              Exclude complete 3'UTR sequence from output.  
-a, --annotation       Exclude gene and GO annotation from output.  
-m, --utrscan_motifs   Exclude UTRscan results from output.  
-t, --teiresias_motifs Exclude Teiresias results from output.  
-s, --target_summary  Exclude summary on target transcripts, genes and GOs from output.  
-r, --mirna            Exclude miRNA target prediction results from output.

### Output Structure:
Header, simply lists the input gene symbol.
```
Input Gene Symbol: mavs
```
Summary data, the first subsection consists of a line with the total number of mRNA transcripts annotated with the query, followed by corresponding lists of Accession numbers, and SeqIDs for each transcript, in the same order. The second subsection shows the number of predicted miRNA response elemnts (MREs) in all mRNA transcripts annotated with the query gene, followed by a non-redundant list of all miRNAs predicted to target one of the mRNAs annotated with the query. The third subsection shows the total number of unique GO terms associated with all the target miRNAs, and corresponding non-redundant list of unique GO terms and GO IDs, in the same order. Summary can be removed from the output using the flag -s.
```
Number of Annotated Transcripts: 6
Transcript Accessions: GIYK01000002;GIYK01000003;...
Transcript SeqIDs in TSA: SS9.1;SS9.2;...

Number of predicted MREs for Gene: 5
miRNAs predicted to bind MREs: ssa-miR-16bc-3p;ssa-miR-22a-3p;...

Number of Unique Annotated GO Terms: 7
GO Terms: activation of innate immune response;mitochondrion;integral component of membrane;...
GO IDs: GO:0002218;GO:0005739;...
```
mRNA header, each individual mRNA entry begins with a header listing the number of the entry, followed by the Accession number, SeqID, and and length of the 3'UTR for the transcript in the entry.
```
Transcript 1:

Accession Number; GIYK01000002
SeqID in TSA: SS9.1
3'UTR length: 226
```
Complete 3'UTR sequence, can be removed from output using the argument -u.
```
3'UTR Sequence: CAACCAGGGTTA...
```
Annotation, lists the gene and transcript annotation given by SQANTI, followed by the the gene symbol and description generated from the Blast analysis performed by Omicsbox, followed by the GO terms and IDs from the OmicsBox annotation. GO terms and IDs are provived as comma-separated lists, ordered with corresponding terms and IDs having the same position in their respective lists. Annotation can be removed from the output using the argument -a.
```
SQANTI Annotation Gene Symbol: LOC100316613
SQANTI Annotation Transcript ID: XM_014131490.1

BLAST Annotation Gene Symbol: mavs
BLAST Annotation Description: interferon promoter stimulating protein 1 isoform X1

GO Terms: activation of innate immune response,mitochondrion,integral component of membrane,...
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
AAATGGCAC [28,37] 0.875
AAGTAGTGC [181,190] 0.6428571428571429
```

miRNA target prediction, consists of a line listing the total number of miRNAs predicted to target the mRNA transcript in the current entry, followed by a list of all the miRNAs. For each predicted miRNA response element (MRE), the output then lists the number of the entry, followed by the name of the targeting miRNA and its mature sequence, folloqws by a list of prediction tools supporting the MRE prediction. After this, the output lists the minimum free energy for the miRNA/mRNA binding predicted by RNAhybrid, the position of the MRE in the 3'UTR, and a four line ascii illustration if the predicted miRNA/mRNA binding. The miRNA target prediction can be removed from the output by using the argument -r
```
Number of predicted MREs for SS9.1: 5
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

## GO_ID_Search.py
### Usage

```python GO_ID_Search.py [-h] [-q QUERY] [-i QUERY_FILE] [-p PREFIX] [-c COMPLEXITY] [-u] [-a] [-m] [-t] [-s] [-r] ```  
Alternatively, arguments can be given interactively by running in script and writing the arguments in the provided field as you would for command line.
All outputs are stored in the OUTPUT folder.

### Arguments:
-h, --help Shows help message and exits

Search terms can be provided as an argument, or by giving a file with one search term per line, or both, in which case all terms from both inputs will be used. gene symbols are case sensitive, and must be written as given in other search outputs or in the SQANTI_OmicsBox_Annotation.tsv file in the DATA folder:  
-q QUERY, --query QUERY                Query IDs, one or more GO IDs separated by semicolons, e.g -q GO:0000050;GO:0004053  
-i QUERY_FILE, --query_file QUERY_FILE Query filename. Must be a list of GO IDs, one per line. At least one of Input_File or Query is required.  

Misc parameters:  
-p PREFIX, --prefix                    PREFIX Optional output prefix. If a prefix is provided, outputs will be stored using the name prefix.query.txt, otherwise they will simply be stored as query.txt.  
-c COMPLEXITY, --complexity COMPLEXITY A demical number indicating minimum Trifonov Linguistic Complexity for Teiresias motifs to be included in the output. If not provided, the program defaults to 0.27.  

Output filters. The following arguments can be provided to selectively remove sections of the output file, as described in the section below.  
-u, --utr              Exclude complete 3'UTR sequence from output.  
-a, --annotation       Exclude gene and GO annotation from output.  
-m, --utrscan_motifs   Exclude UTRscan results from output.  
-t, --teiresias_motifs Exclude Teiresias results from output.  
-s, --target_summary  Exclude summary on target transcripts, genes and GOs from output.  
-r, --mirna            Exclude miRNA target prediction results from output.  

### Output Structure:
Header, lists the input GO ID and the associated GO term.
```
Input GO ID: GO:0004842
GO Term: ubiquitin-protein transferase activity
```
Summary data, the first subsection consists of two lines with the total number of mRNA transcripts annotated with the query GO, and the number of unique different genes annotated to those mRNAs, followed by corresponding lists of Accession numbers, SeqIDs, and gene annotation for each transcript, in the same order. Gene annotation defaults to the OmicsBox annotation if present, followed by the SQANTI annotation, or lists --NA-- if no annotation was found. The second subsection shows the number of predicted miRNA response elemnts (MREs) in all mRNA transcripts annotated with the query GO ID, followed by a non-redundant list of all miRNAs predicted to target one of the mRNAs annotated with the query. Summary can be removed from the output using the flag -s.
```
Number of Annotated Transcripts: 421
Number of unique target genes: 136
Transcript SeqIDs in TSA: CG132.2;CG3154.2;...
Transcript SeqIDs in TSA: SS9.1;SS9.2;...
Transcript Annotated Genes: --NA--;A0A3B4EWV1;...

Number of predicted MREs for GO: 579
miRNAs predicted to bind MREs: ssa-let-7a-1-3-3p;ssa-let-7a-2-3p;...
```
mRNA header, each individual mRNA entry begins with a header listing the number of the entry, followed by the Accession number, SeqID, and and length of the 3'UTR for the transcript in the entry.
```
Transcript 1:

Accession Number; GIYK01066742
SeqID in TSA: CG132.2
3'UTR length: 782
```
Complete 3'UTR sequence, can be removed from output using the argument -u.
```
3'UTR Sequence: ACCTACCCACA...
```
Annotation, lists the gene and transcript annotation given by SQANTI, followed by the the gene symbol and description generated from the Blast analysis performed by Omicsbox, followed by the GO terms and IDs from the OmicsBox annotation. GO terms and IDs are provived as comma-separated lists, ordered with corresponding terms and IDs having the same position in their respective lists. Annotation can be removed from the output using the argument -a.
```
SQANTI Annotation Gene Symbol: LOC100316613
SQANTI Annotation Transcript ID: XM_014131490.1

BLAST Annotation Gene Symbol: mavs
BLAST Annotation Description: interferon promoter stimulating protein 1 isoform X1

GO Terms: activation of innate immune response,mitochondrion,integral component of membrane,...
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
AAATACATTT [346,356] 0.5
AATACATTTT [347,357] 0.4375
```

miRNA target prediction, consists of a line listing the total number of miRNAs predicted to target the mRNA transcript in the current entry, followed by a list of all the miRNAs. For each predicted miRNA response element (MRE), the output then lists the number of the entry, followed by the name of the targeting miRNA and its mature sequence, folloqws by a list of prediction tools supporting the MRE prediction. After this, the output lists the minimum free energy for the miRNA/mRNA binding predicted by RNAhybrid, the position of the MRE in the 3'UTR, and a four line ascii illustration if the predicted miRNA/mRNA binding. The miRNA target prediction can be removed from the output by using the argument -r
```
Number of predicted MREs for SS9.1: 5
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

## DATA file structure
### SQANTI_OmicsBox_Annotation.tsv
A Tab Separated Value file containing the annotation of the transcripts in the FL-mRNA transcriptome. Descripton of the entries:  
Accession_number: Accession number of the transcript in TSA.  
SeqID: A unique identifier that provides information on the clustering process for the transcripts. The two letter prefix indicates the method of clustering, and the genome used for SQANTI annotation. SS = cDNA_cupcake with the Salmo Salar genome, ST = cDNA_cupcake with the Salmo Trutta genome, CG = Cogent. The first number indicates the locus the transcript was mapped to. Transcripts with the same prefix and locus numbers have a high likelihood of being spice variants of the same gene. The second number is a unique identifier within each locus.  
Sqanti_Gene: The gene symbol assigned to the transcript by SQANTI, using the genome indicated in the SeqID prefix.  
Sqanti_Transcript: The accesion number assigned to the transcrip by SQANTI, using the genome indicated in the SeqID prefix.  
Gene_OmicsBox: The gene symbol assigned to the transcript by Blast analysis performed by OmicsBox.  
Sequence_Description: The description generated by natural language processing of Blast analysis performed by OmicsBox.  
Annotation_GO_ID: List of GO IDs assigned to the transcripts by OmicsBox, separated by commas.  
Annotation_GO_Term: List of GO terms assigned to the transcripts by OmicsBox, in the same order as the GO IDs, separated by commas.  
FLNCs: The number of separate PacBio sequencing reads supporting the mRNA transcript.  
Sample_IDs: The IDs of the tissues and conditions the transcripts were observed in.  
```
Accession_number	SeqID	Sqanti_Gene	Sqanti_Transcript  Gene_OmicsBox	Sequence_Description	Annotation_GO_ID	Annotation_GO_Term	FLNCs	Sample_IDs
GIYK01000001	SS4.2		novel	ATP6	ATP synthase F0 subunit 6	GO:0005743,GO:0015078,GO:0015986,GO:0016021,GO:0045263	mitochondrial inner membrane,proton transmembrane transporter activity,ATP synthesis coupled proton transport,integral component of membrane,proton-transporting ATP synthase complex coupling factor F(o)	51223	GiU1,GiU4,GiU7,HKU1,HKU4,HKU7,LiU1,LiU4,LiU7,SAV_Challenge,SAV_Control
GIYK01000002	SS9.1	LOC100316613	XM_014131490.1	mavs	interferon promoter stimulating protein 1 isoform X1	GO:0002218,GO:0005739,GO:0016021,GO:0032481,GO:0048525,GO:0051092,GO:0051607	activation of innate immune response,mitochondrion,integral component of membrane,positive regulation of type I interferon production,negative regulation of viral process,positive regulation of NF-kappaB transcription factor activity,defense response to virus	54	GiU1,GiU4,GiU7,HKU1,HKU4,HKU7,SAV_Control
```
### miRNAome.fa
The full sequence of all the mature miRNAs used for target prediction analysis, in fasta format, with the miRNA name in the header.
```
>ssa-let-7a-1-3-3p
CTGTACAGCCTCCTAGCTTTCC
>ssa-let-7a-2-3p
CTATACAACTTACTGTCTTTCC
...
```
### mRNA_3UTR.fasta
The complete 3'UTR sequence for all transcripts in the FL-mRNA transcriptome in TSA submission GIYK01000000, based on CDS prediction performed using TransDecoder, in fasta format. The header lists the Accession number and the SeqID of the sequence, separated by a semi-colon.
```
>GIYK01000001;SS4.2
TGGCACACCAAGCACACGCATACCACATGGTTGA...
```
### CD-Hit_Clusters.txt
Lists the clusters of 3'UTR sequences generated by CD-Hit to reduce analysis load, one cluster per line, with the SeqIDs of sequences similar enough to be clustered together separated by semi-colons, and the first SeqID in each line being used as the representative for that cluster. If a line contains only one entry, the 3'UTR of that transcript was considered unique in the FL-mRNA transcriptome.
```
SS4.2
SS9.2;SS9.1;SS9.4;SS9.5;SS9.6;SS9.7
```
### RNAhybrid_target_prediction_part_X
A comma separated value file containing the RNAhybrid target prediction results, with each MRE prediction as a single line entry separated by colons. The colums indicate, in order: The input 3'UTR cluster, the length of the 3'UTR, the miRNA predicted to target the MRE, the length of the mature miRNA, the minimum free energy calculated for the mRNA/miRNA interaction, the p-value for the MRE (calculated based on the assumption that the sequences are human, not accurate), the position of the MRE in the 3'UTR, and a four-line ascii illustration of the mRNA/miRNA interacted, separated into four colums.
```
SS38.5;SS38.18:2841:ssa-let-7a-1-3-3p:22:-29.8:0.004863:1791:U     AGG   CA          U: GGAGG   GGG  GGGCUGUACA : CUUUC   UCC  UCCGACAUGU :C     GA                C
SS38.10;SS38.11:3109:ssa-let-7a-1-3-3p:22:-29.8:0.005623:2057:U     AGG   CA          U: GGAGG   GGG  GGGCUGUACA : CUUUC   UCC  UCCGACAUGU :C     GA                C
```
### RNAhybrid_plus_2.txt
A tab separated value file showing which prediction tools support a given mRNA/miRNA interaction for all MREs supported by RNA-hybrid and at least two other prediction tools. The columns indicate, in order: The miRNA predicted to target the MRE, the mRNA cluster containing the MRE, the number of prediction tools supporting the MRE, and a list of the prediction tools separated by commas and a space.
```
ssa-miR-218b-5p	SS581.13	4	PITA, miRanda, TargetSpy, RNAhybrid
ssa-miR-23c-5p	SS1348.7	4	PITA, miRanda, TargetSpy, RNAhybrid
```
### Teiresias_k1000_prob_cutoff_5.txt
Teiresias prediction results for all predicted motifs with occurence more than five times over what is expected randomly by the nucleotide distribution. Each line cotains the total number of occurences of the motif in the 3'UTRome, followed by a tab, followed by the number of different 3'UTR clusters the motif was found it, followed by a tab, followed by the sequence of the motif, and a series of number pairs indicating the cluster and position of all occurences of the motif, separated by spaces. The first number in a pair indicates the 3'UTR using 0-based conting, based on the order in CD-Hit_Clusters.txt, so 0 = SS4.2, 1 = SS9.2;SS9.1;SS9.4;SS9.5;SS9.6;SS9.7, 2 = SS14.1, etc. The second number in a pair indicates the start position of the motif in the 3'UTR, again using 0 based counting, so 0 is the first base, 1 is the second, etc.
```
3350	2047	AGTGCACTA 1 184 11 64 12 167 14 400...
```
### uscan_output.txt
UTRScan prediction results. The first 52 lines contain a header describing all the motif classes. The remaining lines all contain an entry for a single predicted motif, containing the 3'UTR cluster it was found in, the length of the representative sequence, the class of the motif, the motifs start and end position in brackets, and the entire motif sequence
```
CG193.3 1893 : 15-LOX-DICE [101,134] : CTCTACCCCC CACT ACG  CCCTTGCTCT GAGC AGG
CG1061.2 1441 : 15-LOX-DICE [776,809] : CCCCCCCTAC TCTC AGG  CCCCCCCTAC TCTC AGG
```
