#!/usr/bin/env python

"""

                                    User Commands

\033[1mSYNPOPSIS\033[0m:
        Parses a .VCF file to extract the genotype of all individuals based on a 
        per-individual coverage threshold. 
        
        It produces another .VCF file with the locus that passed the per-individual
        coverage threshold and it also produces a text file with the genotype of 
        each individual for these loci.

\033[1mUSAGE\033[0m:
        %program <in_vcf> <per-sample_cov_threshold> <num_samples> <samp_compromise> 
        <output_parsed_vcf> <output_genotype_file>
        
        samp_compromise = number of samples that will be taken off the total number
        of samples in which the coverage threshold has to be respected.
        
\033[1mOUTPUT GENOTYPE FORMAT\033[0m:
        ------------------------------------------------------------------------------
        Contig        Pos    Pop1    Pop2     Ind1    Ind2    Ind3    Ind4    ...    
        contig_12     356    8/16    24/0     A/T     A/A      A/T     T/T    ...    
        ...
        
        NB : The popualtion counts are in numbers of alleles, considering DIPLOID samples.
        It means that if there is 12 samples/pop, the total count per population, if all 
        the samples have been successfully sequenced, is 24.

\033[1mDESCRIPTION OF THE FILTERS\033[0m:
        In order to keep a SNP and output it in the .VCF file, it has to pass the 
        coverage and quality filters, i.e it has to be present in at least <num_samples>
        minus <samp_compromise> individuals with a per-individual coverage chosen by the 
        user (<cov_threshold>). Plus, SNPs with a quality of <= 20 are discarded.

\033[1mCREDITS\033[0m:
        Doctor Pants 2012 \m/

"""

import sys

try:
    in_vcf = sys.argv[1]
    cov_threshold = int(sys.argv[2])
    num_samples = int(sys.argv[3])
    samp_compromise = int(sys.argv[4])
    output_vcf = sys.argv[5]
    output_genotypes = sys.argv[6]
except:
    print __doc__
    sys.exit(2)

# Parses the input VCF file and outputs the loci that pass the per-individual cov 
# threshold in another VCF file specified in the command line (<output_parsed_vcf>).

samples_good_cov = set()

with open(in_vcf, "rU") as vcf_i:
    with open(output_vcf, "w") as vcf_o:
        for line in vcf_i: 
            
            # Copies the vcf header in the vcf parsed output
            if line.startswith("#"):
                line = line.strip()
                vcf_o.write(line + "\n")
            
            if line.startswith("k24"):
                line = line.strip()
                
                # Discards the loci with a SNP quality < 20
                if float(line.split("\t")[5]) < 20:
                    pass
                
                # Discards the major INDELs but keeps the micro-INDELs
                elif line.split()[7].split(";")[0] == "INDEL":
                    
                    if len(line.split()[3]) or len(line.split()[4]) > 2:
                        pass
                    
                    elif len(line.split()[3]) and len(line.split()[4]) <= 2:
                        
                        if int(line.split()[7].split(";")[1].split("=")[1]) >= (cov_threshold*(num_samples-samp_compromise)):
                        
                            sample_count = num_samples
                            ind_column = 9
                            smp_good_cov = 0
                            
                            # Verifies if every sample respects the coverage threshold
                            while sample_count > 0:
                                if int(line.split("\t")[ind_column].split(":")[2]) >= cov_threshold:
                                    smp_good_cov += 1
                                ind_column += 1
                                sample_count -= 1
                            
                            # If there is enough samples that respects the cov. threshold the locus
                            # is kept and written down in the output file.
                            if smp_good_cov == (num_samples - samp_compromise):
                                vcf_o.write(line + "\n")
                        
                # Discards complex SNPs (tri-nucleotide)
                elif len(line.split("\t")[3]) and len(line.split("\t")[4]) > 2:
                    pass
                
                # Keeps the loci with a depth of coverage (DP) >= to the coverage threshold
                # chosen by the user times the number of samples that have to have that 
                # coverage.
                elif int(line.split("\t")[7].split(";")[0].split("=")[1]) >= (cov_threshold*(num_samples-samp_compromise)):
                    
                    sample_count = num_samples
                    ind_column = 9
                    smp_good_cov = 0
                    
                    # Verifies if every sample respects the coverage threshold
                    while sample_count > 0:
                        if int(line.split("\t")[ind_column].split(":")[2]) >= cov_threshold:
                            smp_good_cov += 1
                        ind_column += 1
                        sample_count -= 1
                    
                    # If there is enough samples that respects the cov. threshold the locus
                    # is kept and written down in the output file.
                    if smp_good_cov >= (num_samples - samp_compromise):
                        vcf_o.write(line + "\n")
                        
    print "\nParsed vcf file created\nNow genotyping..."

loci = {}
smp_col = {}

# Extracts the genotype of the samples in the vcf parsed file
with open(output_vcf, "rU") as parsed_vcf:
        
    for line in parsed_vcf:
        
        if line.startswith("##"):
            pass
        
        # Associates the name of the sample to its corresponding column position in a dictionnary
        elif line.startswith("#CHROM"):
            line = line.strip()
            samples = num_samples
            ind_column = 9
            col_num = 1
            
            while samples > 0:
                col = "col" + str("%02i" % col_num)
                sp = "C-" + line.split()[ind_column].split(".")[2]
                smp_col[col] = sp
                ind_column += 1
                samples -= 1
                col_num += 1
        
        elif line.startswith("k24"):
            # Creates the dictionnary containing all the information for every contig
            # (SNP position, REF allele, ALT allele, genotype per individual, allele counts per pop.)
            contig = line.split()[0].split("_")[2] + "_" + str("%05i" % int(line.split()[0].split("_")[3]))
            pos = int(line.split()[1])
            loci[contig] = {}
            loci[contig][pos] = {}
            loci[contig][pos]["alleles"] = {}
            loci[contig][pos]["genotypes"] = {}
            loci[contig][pos]["normals"] = {}
            loci[contig][pos]["normals"]["REF"] = 0
            loci[contig][pos]["normals"]["ALT"] = 0
            loci[contig][pos]["num_normals"] = 0
            loci[contig][pos]["dwarfs"] = {}
            loci[contig][pos]["dwarfs"]["REF"] = 0
            loci[contig][pos]["dwarfs"]["ALT"] = 0
            loci[contig][pos]["num_dwarfs"] = 0
            
            # Writes the REF and ALT alleles in the dictionnary for the locus
            loci[contig][pos]["alleles"]["REF"] = line.split()[3]
            loci[contig][pos]["alleles"]["ALT"] = line.split()[4]
            
            samples = num_samples
            ind_column = 9
            col_num = 1
            
            # Goes through every sample column and 
            while samples > 0:
                col = "col" + str("%02i" % col_num)
                # uses the dictionnary containing the sample names to associate the genotype
                # to the good sample according to the column number.
                sp = smp_col[col]
                loci[contig][pos]["genotypes"][sp] = []
                
                # If the sample respects the coverage threshold, the genotype is written in the
                # dictionnary
                if int(line.split()[ind_column].split(":")[2]) >= cov_threshold:
                    
                    # If the cov. threshold is respected, it adds 1 to the number of individuals
                    # that passed the threshold. In other words, it counts the number of samples
                    # that have a sufficient coverage to be genotyped.
                    if sp.startswith("C-N"):
                        loci[contig][pos]["num_normals"] += 1
                    if sp.startswith("C-D"):
                        loci[contig][pos]["num_dwarfs"] += 1
                    
                    # First allele is REF
                    if line.split()[ind_column].split(":")[0].split("/")[0] == "0":
                        
                        loci[contig][pos]["genotypes"][sp].append(loci[contig][pos]["alleles"]["REF"])
                        if sp.startswith("C-N"):
                            loci[contig][pos]["normals"]["REF"] += 1
                        if sp.startswith("C-D"):
                            loci[contig][pos]["dwarfs"]["REF"] += 1
                    
                    # Second allele is ALT
                    if line.split()[ind_column].split(":")[0].split("/")[0] == "1":
                        
                        loci[contig][pos]["genotypes"][sp].append(loci[contig][pos]["alleles"]["ALT"])
                        if sp.startswith("C-N"):
                            loci[contig][pos]["normals"]["ALT"] += 1
                        if sp.startswith("C-D"):
                            loci[contig][pos]["dwarfs"]["ALT"] += 1
                    
                    # First allele is REF
                    if line.split()[ind_column].split(":")[0].split("/")[1] == "0":
                        
                        loci[contig][pos]["genotypes"][sp].append(loci[contig][pos]["alleles"]["REF"])
                        if sp.startswith("C-N"):
                            loci[contig][pos]["normals"]["REF"] += 1
                        if sp.startswith("C-D"):
                            loci[contig][pos]["dwarfs"]["REF"] += 1
                    
                    # Second allele is ALT
                    if line.split()[ind_column].split(":")[0].split("/")[1] == "1":
                        
                        loci[contig][pos]["genotypes"][sp].append(loci[contig][pos]["alleles"]["ALT"])
                        if sp.startswith("C-N"):
                            loci[contig][pos]["normals"]["ALT"] += 1
                        if sp.startswith("C-D"):
                            loci[contig][pos]["dwarfs"]["ALT"] += 1
                
                # If the sample coverage is below threshold, the genotype is "0"
                if int(line.split()[ind_column].split(":")[2]) < cov_threshold:
                    
                    loci[contig][pos]["genotypes"][sp].append("0")
                
                ind_column += 1
                samples -= 1
                col_num += 1

    print "\nGenotypes extracted\nCreating output file..."

# Writes down the infprmation for each contig in the genotype output file
with open(output_genotypes, "w") as out_geno:
    
    # Creates the first part of the file's header
    out_geno.write("Contig" + "\t" + "POS" + "\t" + "REF" + "/" + "ALT" + "\t" + "Dwarfs" + "\t" + "Normals")
    
    # Outputs the samples IN ORDER in the file's header
    order = []
    for col in smp_col:
        order.append(smp_col[col])
    
    for sample in sorted(order):
        out_geno.write("\t" + sample)
    
    # Goes through every contig in the mega-dictionnary and outputs the information per sample and per pop.
    for contig in list(sorted(loci)):
        for pos in list(sorted(loci[contig])):
            
            loci[contig][pos]["normals"]["REF"]
            
            num_a1_D = loci[contig][pos]["dwarfs"]["REF"]
            num_a2_D = loci[contig][pos]["dwarfs"]["ALT"]
            num_a1_N = loci[contig][pos]["normals"]["REF"]
            num_a2_N = loci[contig][pos]["normals"]["ALT"]
            
            num_dwarfs = loci[contig][pos]["num_dwarfs"]
            num_normals = loci[contig][pos]["num_normals"]
            
            #If an allele is very rare, it's discarded
            if num_a1_D < 3 and num_a1_N < 3:
                pass
            
            elif num_a2_D < 3 and num_a2_N < 3:
                pass
            
            # If there is less than samp_compromise/2 per population that is successfully genotyped
            # the locus is discarded
            elif num_dwarfs < (num_samples - samp_compromise)/2 or num_normals \
                < (num_samples - samp_compromise)/2:
                pass
            
            else:
                # Outputs the contig and the position on the contig
                out_geno.write("\n" + contig + "\t" + str(pos))
                
                # Outputs the genotypes for that position
                out_geno.write("\t" + loci[contig][pos]["alleles"]["REF"] + "/" + 
                loci[contig][pos]["alleles"]["ALT"])
                
                # Outputs the population information (allele count/pop)
                out_geno.write("\t" + str(loci[contig][pos]["dwarfs"]["REF"]) + "/" + 
                str(loci[contig][pos]["dwarfs"]["ALT"]))
                out_geno.write("\t" + str(loci[contig][pos]["normals"]["REF"]) + "/" + 
                str(loci[contig][pos]["normals"]["ALT"]))
                
                # Outputs the genotype for each sample
                for sample in list(sorted(loci[contig][pos]["genotypes"])):
                        
                    if loci[contig][pos]["genotypes"][sample][0] != "0":
                        out_geno.write("\t" + loci[contig][pos]["genotypes"][sample][0] + "/" + 
                        loci[contig][pos]["genotypes"][sample][1])
                        
                    if loci[contig][pos]["genotypes"][sample][0] == "0":
                        out_geno.write("\t" + loci[contig][pos]["genotypes"][sample][0])

    print "Genotype file successfully created\n"
