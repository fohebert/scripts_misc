#!/usr/bin/python

"""
v.0.1                           User's Commands                         v.0.1

\033[1mDESCRIPTION\033[0m
    Input files: 
    1) edgeR/limma topTags table; 
    2) clusters.txt (Corset output); 
    3) trinotate report
    -------------------------------------------------
    Output files:
    1) TSV file with Transcript-Cluster-GeneID
    2) edgeR/limma topTags table with GeneID added

    This script merges outputs from edgeR/limma, trinotate
    and Corset by adding a column to the topTags file right
    beside the first column so that one knows exactly to
    what protein/gene each sequence corresponds. So when a
    sequence is diff. expressed between 2 conditions, it's
    easy to see what this gene is (e.g. phosphatase). It
    also produces a simple TSV file that contains, for
    each Trinity assembled contig, the name of the corset
    cluster to which it correspond and all the annotation
    terms found for this cluster.

\033[1mUSAGE\033[0m
    %program <DGE_table> <clusters> <trinotate> <out_DGE_table> <out_TSV_annotation>

\033[1mCREDITS\033[0m
    Doc Pants 2015 \m/
"""

import sys
import re

try:
    DGE_table = sys.argv[1]
    corset_clusters = sys.argv[2]
    trinotate_report = sys.argv[3]
    out_DGE_table = sys.argv[4]
    out_TSV_annotation = sys.argv[5]
except:
    print __doc__
    sys.exit(2)

# GLOBAL VARIABLES
# These will be used throughout the whole script
clusters = {} # contains ALL the info (transcripts, DGE info, annotation) - ***META DICTIONNARY***
DGE_in_header = {} # matches each column num. in DGE file to their corresponding name
trspt_annotation = {} # stores annotation info for each transcript

# Dealing with Trinotate annotation file
with open(trinotate_report, "rU") as trinotate_in:
    for line in trinotate_in:
        line = line.strip()
        
        annotation = ""
        GO_terms = []
        # Fills in the annotation section of the meta dict() (clusters)
        if line.startswith("#gene_id") == False:
            trspt_id = line.split("\t")[1] # transcript name
            if trspt_id not in trspt_annotation:
                
                # Adding transcript name into trspt_annotation dict()
                trspt_annotation[trspt_id] = {} # Builds dict() for each transcript in 'trspt_annotation'

                # *** ANNOTATION ***
                # 1) Top BLASTx annotation comes first
                if line.split("\t")[2] != ".":
                    annotation = line.split("\t")[2].split("^")[0] + "^" + line.split("\t")[2].split("^RecName: Full=")[-1].split(";")[0]
                
                # 2) Top BLASTp comes second
                if line.split("\t")[7] != ".":
                    blastp_anno = line.split("\t")[7].split("^")[0] + "^" + line.split("\t")[7].split("^RecName: Full=")[-1].split(";")[0]
                    if blastp_anno != annotation:
                        annotation += ":" + blastp_anno
                
                # 3) If annotation is still empty, will look at pFAM and signalIP annotation
                if annotation == "":
                    # Pfam annotation
                    if line.split("\t")[9] != ".":
                        annotation += line.split("\t")[9].split("^")[2]
                    # signalIP annotation
                    elif line.split("\t")[10] != ".":
                        if line.split("\t")[10].split("^")[-1] == "YES":
                            annotation += "unknown secreted protein"
                
                # *** GO terms ***
                # Combines all 3 sources of GO terms (eggnog,blast,pfam)
                if line.split("\t")[12] != ".":
                    # EGGnog annotation
                    GO_terms.append(line.split("\t")[12].split("^")[-1])
                if line.split("\t")[13] != ".":
                    # Gene Ontology - BLAST+
                    GO_terms.append(line.split("\t")[13])
                if line.split("\t")[14] != ".":
                    # Gene Ontology - pFAM
                    GO_terms.append(line.split("\t")[14])
                
                # If, after all these steps, the annotation and GO term variables are still empty
                # it will fill them with a dot
                if annotation == "":
                    annotation = "."
                if GO_terms == []:
                    GO_terms.append(".")
                
                # Fills in the dictionnary with corresponding annotation + GO terms (empty or not)
                trspt_annotation[trspt_id]["annotation"] = annotation
                trspt_annotation[trspt_id]["GO"] = GO_terms

# Dealing with Corset cluster file
with open(corset_clusters, "rU") as clst_in:
    # Each cluster name is entered in meta dict() (clusters) with its respective transcripts
    # that form the cluster
    for line in clst_in:
        line = line.strip()
        
        # Isolating cluster and transcript names
        clst = line.split("\t")[1]
        transcript = line.split("\t")[0]
        
        # Stores the names of the transcripts that clustered together
        # Each transcript will then be stored with its corresponding Trinotate annotation
        if clst not in clusters:
            clusters[clst] = {}
            clusters[clst]["transcripts"] = {}
            clusters[clst]["transcripts"][transcript] = {}
            clusters[clst]["transcripts"][transcript]["annotation"] = trspt_annotation[transcript]["annotation"]
            clusters[clst]["transcripts"][transcript]["GO"] = trspt_annotation[transcript]["GO"]
        else clst in clusters:
            clusters[clst]["transcripts"][transcript] = {}
            clusters[clst]["transcripts"][transcript]["annotation"] = trspt_annotation[transcript]["annotation"]
            clusters[clst]["transcripts"][transcript]["GO"] = trspt_annotation[transcript]["GO"]

# Dealing with DGE table from limma or EdgeR
with open(DGE_table, "rU") as in_DGE:
    for line in in_DGE:
        line = line.strip()
        
        # Keeps in a dictionnary the name of each column in DGE table
        if re.findall("P.Value", line):
            tot_num_col = len(line.split("\t"))
            count = 1
            while count <= tot_num_col:
                DGE_in_header[count] = line.split("\t")[count-1]
                count += 1

        # When the program finds any line after the header, it keeps the info
        # in the meta dictionnary
        elif line.startswith("Cluster-"):
            # Gets the name of the cluster in DGE table
            clst = line.split("\t")[0]
            # Fills in meta dict() with DGE info: logFC, P.value, Adj.P.Val/FDR, etc. for each cluster
            tot_num_col = len(line.split("\t")) - 1
            count = 1
            while count <= tot_num_col:
                col_name = DGE_in_header[count]
                clusters[clst][col_name] = line.split("\t")[count]
                count += 1

# Creating output files
with open(out_DGE_table, "w") as out_DGE:
    with open(out_TSV_annotation, "w") as out_anno:
        
        # Creates DGE table header
        out_DGE.write("\t" + "Annotation")
        for col_num in DGE_in_header:
            out_DGE.write("\t" + DGE_in_header[col_num])
        
        # Fills in DGE table, cluster by cluster
        for clst in sorted(clusters):
            
            # Gathers annotation info for all the transcripts in a given cluster
            tmp_anno = [] # Place anno. info in temporary list
            for trspt in clusters[clst]["transcripts"]:
                tmp_anno.append(clusters[clst]["transcripts"][trspt]["annotation"])
                
                # Writing down transcript annotation + GO information in annotation output file
                out_anno.write(trspt + "\t" + clst + "\t" + clusters[clst]["transcripts"][trspt]["annotation"] + "\t" + "`".join(clusters[clst]["transcripts"][trspt]["GO"]) + "\n")

            if len(clusters[clst]) > 1:
                # Writes the name of the cluster in DGE output file
                out_DGE.write("\n" + clst)
            
                # Joins pieces of annotation info together in 1 object
                col_to_add_to_DGE = "^".join(tmp_anno)
            
                # Writes down "merged" annotation info as first column of "new" DGE table
                out_DGE.write("\t" + col_to_add_to_DGE)
            
                # Writes down rest of DGE table in output file (e.g. logFC, P.Value, etc.)
                for col_num in DGE_in_header:
                    out_DGE.write("\t" + clusters[clst][DGE_in_header[col_num]])
            
print "\n\033[1mJob Done\033[0m\n"
