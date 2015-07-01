#!/usr/bin/env python

"""

                                    User Commands

\033[1mSYNPOPSIS\033[0m:
        Parses a .VCF file to extract the genotype of all the samples based on  a
        per-sample coverage threshold. It returns in an output file the genotypes
        at the loci for which a sufficient number of samples were  sequenced with
        a minimum depth of coverage that is chosen by the user.

\033[1mUSAGE\033[0m:
        %program <vcf> <per-sample_cov_threshold> <num_samples> <samp_compromise> 
        <output_genotype_file>
        
        samp_compromise = maximum  number  of  samples allowed to be discarded in
        total in  order to keep the locus in the analysis. This number divided by
        2 is the per-population  threshold,  i.e the number of individuals in one
        population that can be discarded because of insufficient coverage.
        
\033[1mOUTPUT GENOTYPE FORMAT\033[0m:
        --------------------------------------------------------------------------
        Contig        Pos    Pop1    Pop2     Ind1    Ind2    Ind3    Ind4    ...    
        contig_12     356    8/16    24/0     A/T     A/A      A/T     T/T    ...    
        ...
        
        NB : The  popualtion counts are in numbers of alleles, considering DIPLOID 
        samples. It means  that  if  there  is 12 samples/pop, the total count per 
        population, if all the samples have been successfully sequenced, is 24.

\033[1mDESCRIPTION OF THE FILTERS\033[0m:
        In order to keep a SNP and output it in the final file, it  has to pass the 
        coverage and quality filters,  i.e  it  has  to  be  present  in  at  least 
        <num_samples> minus <samp_compromise> individuals with a per-individual co-
        verage chosen by the user  (<cov_threshold>). Plus, SNPs with a quality  of 
        <= 20 are discarded.

\033[1mCREDITS\033[0m:
        Doctor Pants 2012 \m/

"""

import sys

try:
    in_vcf = sys.argv[1]
    cov_threshold = int(sys.argv[2])
    num_samples = int(sys.argv[3])
    samp_compromise = int(sys.argv[4])
    output_genotypes = sys.argv[5]
except:
    print __doc__
    sys.exit(2)

loci = {}
smp_col = {}
coverage = {}
complex_SNP = 0
bad_quality = 0
too_rare_allele = 0
not_enough_coverage = 0
SNPs = 0
# Extracts the genotype of the samples in the vcf parsed file
with open(in_vcf, "rU") as vcf_i:
    
    for line in vcf_i:
        
        if line.startswith("##"):
            pass
        
        # Associates the name of the sample to its corresponding column position in a dictionnary
        elif line.startswith("#CHROM"):
            line = line.strip()
            samples = num_samples
            ind_column = 9
            
            while samples > 0:
                sp = "C-" + line.split("\t")[ind_column].split(".")[3]
                smp_col[ind_column] = sp
                coverage[sp] = []
                
                ind_column += 1
                samples -= 1
        
        elif line.startswith("#") == False:
        
            # Creates the dictionnary containing all the pieces of information for every gene
            # (SNP position, REF allele, ALT allele, genotype per individual, allele counts per pop.)
            gene = line.split("\t")[0]
            pos = int(line.split("\t")[1])
            SNPs += 1
            
            if len(line.split("\t")[3]) and len(line.split("\t")[4]) > 1:
                    complex_SNP += 1
                    pass
            
            elif float(line.split("\t")[5]) < float(20):
                    bad_quality += 1
                    pass
            
            else:
            
                if gene not in loci:
                    loci[gene] = {}
                
                loci[gene][pos] = {}
                loci[gene][pos]["alleles"] = {}
                loci[gene][pos]["genotypes"] = {}
                loci[gene][pos]["normals"] = {}
                loci[gene][pos]["normals"]["coverage"] = 0
                loci[gene][pos]["normals"]["REF"] = 0
                loci[gene][pos]["normals"]["ALT"] = 0
                loci[gene][pos]["num_normals"] = 0
                loci[gene][pos]["dwarfs"] = {}
                loci[gene][pos]["dwarfs"]["coverage"] = 0
                loci[gene][pos]["dwarfs"]["REF"] = 0
                loci[gene][pos]["dwarfs"]["ALT"] = 0
                loci[gene][pos]["num_dwarfs"] = 0
                
                # Writes the REF and ALT alleles in the dictionnary for the locus
                loci[gene][pos]["alleles"]["REF"] = line.split("\t")[3]
                loci[gene][pos]["alleles"]["ALT"] = line.split("\t")[4]
                
                samples = num_samples
                ind_column = 9
                
                # Goes through every sample column and verifies if the cov. threshold is respected
                # and if it is, the genotype is extracted for this sample.
                while samples > 0:
                    # uses the dictionnary containing the sample names to associate the genotype
                    # to the good sample according to the column number.
                    sp = smp_col[ind_column]
                    loci[gene][pos]["genotypes"][sp] = []
                    
                    # Fills in the dictionnary "coverage" that contains the information on the coverage
                    # for each sample at each locus. Will be used later to calculate the mean coverage
                    # per sample.
                    cov = int(line.split("\t")[ind_column].split(":")[-3])
                    coverage[sp].append(cov)
                    
                    # If the sample respects the coverage threshold, the genotype is written in the
                    # dictionnary
                    if int(line.split("\t")[ind_column].split(":")[-3]) >= cov_threshold:
                        
                        # If the cov. threshold is respected, it adds 1 to the number of individuals
                        # that passed the threshold. In other words, it counts the number of samples
                        # that have a sufficient coverage to be genotyped.
                        if sp.startswith("C-N"):
                            loci[gene][pos]["num_normals"] += 1
                            loci[gene][pos]["normals"]["coverage"] += int(line.split("\t")[ind_column].split(":")[-3])
                        if sp.startswith("C-D"):
                            loci[gene][pos]["num_dwarfs"] += 1
                            loci[gene][pos]["dwarfs"]["coverage"] += int(line.split("\t")[ind_column].split(":")[-3])
                        
                        # First allele is REF
                        if line.split("\t")[ind_column].split(":")[0].split("/")[0] == "0":
                            
                            loci[gene][pos]["genotypes"][sp].append(loci[gene][pos]["alleles"]["REF"])
                            if sp.startswith("C-N"):
                                loci[gene][pos]["normals"]["REF"] += 1
                            if sp.startswith("C-D"):
                                loci[gene][pos]["dwarfs"]["REF"] += 1
                        
                        # Second allele is ALT
                        if line.split("\t")[ind_column].split(":")[0].split("/")[0] == "1":
                            
                            loci[gene][pos]["genotypes"][sp].append(loci[gene][pos]["alleles"]["ALT"])
                            if sp.startswith("C-N"):
                                loci[gene][pos]["normals"]["ALT"] += 1
                            if sp.startswith("C-D"):
                                loci[gene][pos]["dwarfs"]["ALT"] += 1
                        
                        # First allele is REF
                        if line.split("\t")[ind_column].split(":")[0].split("/")[1] == "0":
                            
                            loci[gene][pos]["genotypes"][sp].append(loci[gene][pos]["alleles"]["REF"])
                            if sp.startswith("C-N"):
                                loci[gene][pos]["normals"]["REF"] += 1
                            if sp.startswith("C-D"):
                                loci[gene][pos]["dwarfs"]["REF"] += 1
                        
                        # Second allele is ALT
                        if line.split()[ind_column].split(":")[0].split("/")[1] == "1":
                            
                            loci[gene][pos]["genotypes"][sp].append(loci[gene][pos]["alleles"]["ALT"])
                            if sp.startswith("C-N"):
                                loci[gene][pos]["normals"]["ALT"] += 1
                            if sp.startswith("C-D"):
                                loci[gene][pos]["dwarfs"]["ALT"] += 1
                    
                    # If the sample coverage is below threshold, the genotype is "0"
                    if int(line.split("\t")[ind_column].split(":")[-3]) < cov_threshold:
                        
                        loci[gene][pos]["genotypes"][sp].append("0")
                    
                    ind_column += 1
                    samples -= 1

    print "\nGenotypes extracted\nCreating output file..."

# Writes down the information for each contig in the genotype output file
with open(output_genotypes, "w") as out_geno:
    
    # Creates the first part of the file's header
    out_geno.write("Gene" + "\t" + "POS" + "\t" + "REF" + "/" + "ALT" + "\t" + "cov_dwarf" + "\t" + "cov_normal" + "\t" + "Dwarfs" + "\t" + "Normals")
    
    # Outputs the samples IN ORDER in the file's header
    order = []
    for col in smp_col:
        order.append(smp_col[col])
    
    for sample in sorted(order):
        out_geno.write("\t" + sample)
    
    # Goes through every contig in the mega-dictionnary and outputs the information per sample and per pop.
    for gene in list(sorted(loci)):
        
        for pos in list(sorted(loci[gene])):
            
            num_a1_D = loci[gene][pos]["dwarfs"]["REF"]
            num_a2_D = loci[gene][pos]["dwarfs"]["ALT"]
            num_a1_N = loci[gene][pos]["normals"]["REF"]
            num_a2_N = loci[gene][pos]["normals"]["ALT"]
            
            num_dwarfs = loci[gene][pos]["num_dwarfs"]
            num_normals = loci[gene][pos]["num_normals"]
            
            # If an allele has a minor allele count < 2 across both populations, the locus is discarded
            if num_a1_D + num_a1_N < 2:
                too_rare_allele += 1
                pass
            
            elif num_a2_D + num_a2_N < 2:
                too_rare_allele += 1
                pass
            
            # If there is less than samp_compromise/2 per population that is successfully genotyped
            # the locus is discarded
            elif num_dwarfs < (num_samples - samp_compromise)/2 or num_normals \
                < (num_samples - samp_compromise)/2:
                not_enough_coverage += 1
                pass
            
            else:
                # Outputs the contig and the position on the contig
                out_geno.write("\n" + gene + "\t" + str(pos))
                
                # Outputs the genotypes for that position
                out_geno.write("\t" + loci[gene][pos]["alleles"]["REF"] + "/" + \
                    loci[gene][pos]["alleles"]["ALT"])
                
                out_geno.write("\t" + str(loci[gene][pos]["dwarfs"]["coverage"]))
                out_geno.write("\t" + str(loci[gene][pos]["normals"]["coverage"]))
                
                # Outputs the population information (allele count/pop)
                out_geno.write("\t" + str(loci[gene][pos]["dwarfs"]["REF"]) + "/" + \
                    str(loci[gene][pos]["dwarfs"]["ALT"]))
                out_geno.write("\t" + str(loci[gene][pos]["normals"]["REF"]) + "/" + \
                    str(loci[gene][pos]["normals"]["ALT"]))
                
                # Outputs the genotype for each sample
                for sample in list(sorted(loci[gene][pos]["genotypes"])):
                        
                    if loci[gene][pos]["genotypes"][sample][0] != "0":
                        out_geno.write("\t" + loci[gene][pos]["genotypes"][sample][0] + "/" + \
                            loci[gene][pos]["genotypes"][sample][1])
                        
                    if loci[gene][pos]["genotypes"][sample][0] == "0":
                        out_geno.write("\t" + loci[gene][pos]["genotypes"][sample][0])

    for sample in list(sorted(coverage)):
        print sample, "coverage: ", float(sum(coverage[sample])) / float(len(coverage[sample]))
        with open(sample, "w") as cov_file:
            for snp_cov in coverage[sample]:
                cov_file.write(str(snp_cov) + "\n")

    total = complex_SNP + bad_quality + too_rare_allele + not_enough_coverage
    print "\nComplex SNPs discarded = ", complex_SNP, "\nDiscarded for low quality = ", bad_quality, "\nToo rare = ", too_rare_allele, "\nInsufficient coverage = ", not_enough_coverage
    print "\nTotal discarded/total SNPs = ", total,"/",SNPs 
    print "\nGenotype file successfully created\n"
