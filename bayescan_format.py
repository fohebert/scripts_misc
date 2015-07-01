#!/usr/bin/python
# -*- coding: utf-8 -*-

"""

\033[1mSYNOPSIS\033[0m

        Formats the file from the script "SNP_count.py" in bayescan format.
        It takes the genotype information from the CLC SNP table and counts
        the number of allele for each sample. It writes that info in an out-
        put file that respects the input format file for bayescan.
        
        It is \033[1mIMPORTANT\033[0m to use the script:
        
        --> bayescan_loci_to_discard.py <--
        
        Bayescan cannot treat the loci for which a subset of samples have
        been genotyped. It needs only the loci for which the genotype info
        is complete for every sample (if there is 5 samples per pop. and 
        for one locus, only 4 samples are successfully genotyped, it is
        important to discard that locus in the bayescan analysis).

\033[1mUSAGE\033[0m

        %program <input_file> <genotypes> <output_file>
        
        \033[1mIMPORTANT : CHANGE DE NUMBER OF SAMPLES PER POP AND
        THE NUMBER OF LOCI DIRECTLY IN THE SCRIPT !!!!!!!\033[0m
        
        output file format:
        
        __________file starts below the line__________
        
        [loci]=508

        [Populations]=2

        [pop]=1
        1  6  2  2 4 
        2  6  2  2 4 
        3  6  2  2 4 
        4  6  2  2 4 
        ...
        
        [pop]=2
        1  6  2  3 3 
        2  6  2  3 3 
        3  6  2  3 3 
        4  6  2  2 4 
        ... 
        
        ___________file ends above the line___________

\033[1mCREDITS\033[0m

        Doctor Pants 2011 \m/

"""

import sys

try:
    input_file = sys.argv[1]
    genotypes = sys.argv[2]
    out_file = sys.argv[3]
except:
    print __doc__
    sys.exit(0)

samples_per_pop = 3
max_loci = 508

# Stores the genotype information at each locus

locus_genotype = {}

with open(genotypes, "r") as in_f:
    for line in in_f:
        line = line.strip()
        if line != "":
            locus_n = "cocl_" + line.split("g")[1].split()[0] + "_" + line.split()[1]
            genotype = line.split()[5] + line.split()[14]
            locus_genotype[locus_n] = [genotype[0],genotype[1]]

# Constructs the actual output file

dwarf = {}
normal = {}
first_dwarf = False
locus_num = 0
allele1_D = 0
allele2_D = 0
allele1_N = 0
allele2_N = 0

with open(input_file, "r") as in_f:
    for line in in_f:
        line = line.strip()
        if line != "":
                
            # Treast the first line
            if locus_num == 0:
                
                locus_num += 1
                locus_name = line.split()[0]
                A1 = locus_genotype[locus_name][0]
                A2 = locus_genotype[locus_name][1]
                
                if line.split()[2].split("-")[1][0] == "D":
                    
                    first_dwarf = True
                    sample_allele1_D = 0
                    sample_allele2_D = 0
                    
                    # Allele count for each individual: if there is
                    # a homozygous, the count is 0 for one allele 
                    # and 2 for the other.
                    if float(line.split("\t")[3]) > 0:
                        if A1 == "A":
                            sample_allele1_D += 1
                        if A2 == "A":
                            sample_allele2_D += 1
                        
                    if float(line.split("\t")[4]) > 0:
                        if A1 == "C":
                            sample_allele1_D += 1
                        if A2 == "C":
                            sample_allele2_D += 1
                        
                    if float(line.split("\t")[5]) > 0:
                        if A1 == "G":
                            sample_allele1_D += 1
                        if A2 == "G":
                            sample_allele2_D += 1
                        
                    if float(line.split("\t")[6]) > 0:
                        if A1 == "T":
                            sample_allele1_D += 1
                        if A2 == "T":
                            sample_allele2_D += 1
                    
                    if sample_allele1_D == 0:
                        sample_allele2_D += 1
                    if sample_allele2_D == 0:
                        sample_allele1_D += 1
                    
                    # adds the allele count for the individual to
                    # the total allele count for the population at
                    # that locus.
                    allele1_D += sample_allele1_D
                    allele2_D += sample_allele2_D
                        
            # Same thing, but for the other lines.
            if locus_num > 0:
                    
                locus_name = line.split()[0]
                A1 = locus_genotype[locus_name][0]
                A2 = locus_genotype[locus_name][1]   
                                
                # Writes down the info for the locus in treatment in 
                # the input file. 
                if line.split()[2] == "XX_noTag":
                   
                    if 1 <= locus_num < 10:
                    
                        l_num = "00" + str(locus_num)
                    
                        dwarf[int(l_num)] = [(allele1_D), (allele2_D)]
                        normal[int(l_num)] = [(allele1_N), (allele2_N)]
                        
                        locus_num += 1
                    
                    elif 10 <= locus_num < 100:
                    
                        l_num = "0" + str(locus_num)
                    
                        dwarf[int(l_num)] = [(allele1_D), (allele2_D)]
                        normal[int(l_num)] = [(allele1_N), (allele2_N)]
                        
                        locus_num += 1
                    
                    elif locus_num >= 100:
                    
                        dwarf[locus_num] = [(allele1_D), (allele2_D)]
                        normal[locus_num] = [(allele1_N), (allele2_N)]
                        
                        locus_num += 1
                    
                    allele1_D = 0
                    allele2_D = 0
                    allele1_N = 0
                    allele2_N = 0
                
                else:
                    
                    # Treats the other lines here, just like the first one.
                    # This is for the dwarfs
                    if line.split()[2].split("-")[1][0] == "D":
                        
                        if first_dwarf == True:
                            first_dwarf = False
                            continue
                        else:
                            
                            sample_allele1_D = 0
                            sample_allele2_D = 0
                            
                            if float(line.split("\t")[3]) > 0:
                                
                                if A1 == "A":
                                    sample_allele1_D += 1
                                if A2 == "A":
                                    sample_allele2_D += 1
                                
                            if float(line.split("\t")[4]) > 0:
                                
                                if A1 == "C":
                                    sample_allele1_D += 1
                                if A2 == "C":
                                    sample_allele2_D+= 1
                                
                            if float(line.split("\t")[5]) > 0:
                                
                                if A1 == "G":
                                    sample_allele1_D += 1
                                if A2 == "G":
                                    sample_allele2_D += 1
                               
                            if float(line.split("\t")[6]) > 0:
                                
                                if A1 == "T":
                                    sample_allele1_D += 1
                                if A2 == "T":
                                    sample_allele2_D += 1
                            
                            if sample_allele1_D == 0:
                                sample_allele2_D += 1
                            if sample_allele2_D == 0:
                                sample_allele1_D += 1
                            
                            allele1_D += sample_allele1_D
                            allele2_D += sample_allele2_D
                                
                    # This is foro the normals
                    if line.split()[2].split("-")[1][0] == "N":
                        
                        sample_allele1_N = 0
                        sample_allele2_N = 0
                        
                        if float(line.split("\t")[3]) > 0:
                            if A1 == "A":
                                sample_allele1_N += 1
                            if A2 == "A":
                                sample_allele2_N += 1
                            
                        if float(line.split("\t")[4]) > 0:
                            if A1 == "C":
                                sample_allele1_N += 1
                            if A2 == "C":
                                sample_allele2_N += 1
                            
                        if float(line.split("\t")[5]) > 0:
                            if A1 == "G":
                                sample_allele1_N += 1
                            if A2 == "G":
                                sample_allele2_N += 1
                            
                        if float(line.split("\t")[6]) > 0:
                            if A1 == "T":
                                sample_allele1_N += 1
                            if A2 == "T":
                                sample_allele2_N += 1
                        
                        if sample_allele1_N == 0:
                            sample_allele2_N += 1
                        if sample_allele2_N == 0:
                            sample_allele1_N += 1
                            
                        allele1_N += sample_allele1_N
                        allele2_N += sample_allele2_N

with open(out_file, "w") as out_f:
    
    # Constructs the first lines of the output file
    number_of_loci = str(locus_num - 1)
    out_f.write("[loci]=" + number_of_loci + "\n" + "\n" + "[populations]=2" + "\n" + "\n") 
    
    # Writes downs in the output file the num. of observations at each locus for the dwarfs
    out_f.write("[pop]=1" + "\n")
    
    for loci in list(sorted(dwarf)):
        if 1 <= loci < 10:
            out_f.write("  " + str(loci) + "  " + str(samples_per_pop*2) + "  " + "2" + "  ")
        elif 10 <= loci < 100:
            out_f.write(" " + str(loci) + "  " + str(samples_per_pop*2) + "  " + "2" + "  ")
        elif loci >= 100:
            out_f.write(str(loci) + "  " + str(samples_per_pop*2) + "  " + "2" + "  ")
        
        for counts in list(sorted(dwarf[loci])):
        
            out_f.write(str(counts) + " ")
        
        out_f.write("\n")
    
    # Writes downs in the output file the num. of observations at each locus for the normals
    out_f.write("\n" + "[pop]=2" + "\n")
    
    for loci in list(sorted(normal)):
        
        if 1 <= loci < 10:
            out_f.write("  " + str(loci) + "  " + str(samples_per_pop*2) + "  " + "2" + "  ")
        elif 10 <= loci < 100:
            out_f.write(" " + str(loci) + "  " + str(samples_per_pop*2) + "  " + "2" + "  ")
        elif loci >= 100:
            out_f.write(str(loci) + "  " + str(samples_per_pop*2) + "  " + "2" + "  ")
        
        for counts in list(sorted(normal[loci])):
            
            out_f.write(str(counts) + " ")
        
        out_f.write("\n")
    
    print "\nSelf-deceit is an unbearable mistake...\n"
