#!/usr/bin/env python

import gzip
from tqdm import tqdm
from collections import Counter
import numpy as np
from pysam import VariantFile
import pysam



# Proximal_snps function takes args as, ( Joint_vcf file location, O_file_id, F_file_id, M_file_id, list_denovo_region) parameters names are diff but not pos. 
def Proximal_snps(joint_vcf_loc, o_id, f_id, m_id, dnv_list_list):

    # storing all files names in the row into list.
    column_files_indexs = {}
    file_names = []
    with gzip.open(f'{joint_vcf_loc}/1kGP_high_coverage_Illumina.chr1.filtered.SNV_INDEL_SV_phased_panel.vcf.gz', 'rt') as vcf_file:
        
        for line in vcf_file:
            
            line = line.strip().split('\t')
            if line[0].startswith('#CHROM'): # considering only row in header which has file names.
                file_names.extend(line)
                break
    
    
    # getting the index to all files in the row. 
    for file in tqdm(file_names):
        
        file_index = file_names.index(file)
        column_files_indexs[file] = file_index           


    # As our denovo regions file is sort based on chrs, we read the joint vcf file based the current chr in the loop.
    previous_chrom = 'chr1'
    vcf_in = VariantFile('../../../../1000genomes_SNP/1kGP_high_coverage_Illumina.chr1.filtered.SNV_INDEL_SV_phased_panel.vcf.gz')


    # Storing de novo regions which has nearby heterozygous snp, flank size 350.
    NHsnp_denovo_sites = []


    # here we itterate list of denovo regions.
    for Dregion in tqdm(dnv_list_list):
        
        # 
        chrom = Dregion[0]
        if chrom == 'chrX': continue
        if chrom != previous_chrom:
            vcf_in = VariantFile(f'{joint_vcf_loc}1kGP_high_coverage_Illumina.{chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz')
        
        Dnv_chrm = Dregion[0]
        Dnv_start = int(Dregion[1])
        Dnv_end = int(Dregion[2])
        
        
        # Getting the index of the trio files.
        o = o_id
        f = f_id
        m = m_id


        # Calling the file index to fetch the according sample info (GT, Refn, Alt).
        O = column_files_indexs[o]
        F = column_files_indexs[f]
        M = column_files_indexs[m]

        
        # left flank
        snp_alleles_info = {}
        Region_info = vcf_in.fetch(Dnv_chrm, Dnv_start - 350, Dnv_start+1)
        for line in Region_info: # adding the usefull snp info, complete repeat line, reference, alt, gt's, pos - repeat_start
            line = str(line).strip().split('\t')
            # considere only the substitutions.
            if len(line[3]) != len(line[4]) : continue
            Genotype_info = [Dregion,line[O], line[F], line[M], line[1], line[3], line[4], int(Dnv_start) - int(line[1])]
            
            # Considering only bi-allelic loci.
            current_snp_loci = f'{line[0]}-{line[1]}'
            if current_snp_loci in snp_alleles_info:
                snp_alleles_info[current_snp_loci][1] += 1
            else: snp_alleles_info[current_snp_loci] = [Genotype_info, 1]
            
            # keys and values
            keys = list(snp_alleles_info.keys())
            index_0 = keys[0]
            
            #
            if len(snp_alleles_info) > 1: # ********************** final value in the dictionary
                if int(snp_alleles_info[index_0][-1]) > 1:
                    snp_alleles_info.pop(index_0)
                else:               
                    
                    Genotype_info1 = snp_alleles_info[index_0][0] 
                    # if Offspring is hetero or not.
                    if Genotype_info1[1] in ['0|1', '1|0']:

                        # Checking the inheritence from the parents.
                        if Genotype_info1[2] in ['0|1', '1|0'] and Genotype_info1[3] in ['0|0', '1|1']:
                            NHsnp_denovo_sites.append([*Genotype_info1[0], *Genotype_info1[1:], 'Left_flank'])
                            # print(Genotype_info1[0], Genotype_info1[1:], 'Left_flank', file = out, sep = '\t')

                        elif Genotype_info1[3] in ['0|1', '1|0'] and Genotype_info1[2] in ['0|0', '1|1']:
                            NHsnp_denovo_sites.append([*Genotype_info1[0], *Genotype_info1[1:], 'Left_flank'])
                            # print(*Genotype_info1[0], *Genotype_info1[1:], 'Left_flank', file = out, sep = '\t')
                            
                    snp_alleles_info.pop(index_0)  
            
            
           
        # Right flank
        snp_alleles_info_RF = {}
        Region_info = vcf_in.fetch(Dnv_chrm, Dnv_end, Dnv_end + 350 + 1)
        for line in Region_info:
            
            line = str(line).strip().split('\t')
            # considere only the substitutions.
            if len(line[3]) != len(line[4]) : continue
            Genotype_info = [Dregion,line[O], line[F], line[M], line[1], line[3], line[4], int(line[1]) - int(Dnv_end)]
            
            # Considering only bi-allelic loci.
            current_snp_loci = f'{line[0]}-{line[1]}'
            if current_snp_loci in snp_alleles_info_RF:
                snp_alleles_info_RF[current_snp_loci][1] += 1
            else: snp_alleles_info_RF[current_snp_loci] = [Genotype_info, 1]
            
            # keys and values
            keys = list(snp_alleles_info_RF.keys())
            index_0 = keys[0]
            
            #
            if len(snp_alleles_info_RF) > 1:
                if int(snp_alleles_info_RF[index_0][-1]) > 1:
                    snp_alleles_info_RF.pop(index_0)
                else:               
                    
                    Genotype_info1 = snp_alleles_info_RF[index_0][0] 
                    # if Offspring is hetero or not.
                    if Genotype_info1[1] in ['0|1', '1|0']:

                        # Checking the inheritence from the parents.
                        if Genotype_info1[2] in ['0|1', '1|0'] and Genotype_info1[3] in ['0|0', '1|1']:
                            NHsnp_denovo_sites.append([*Genotype_info1[0], *Genotype_info1[1:], 'Right_flank'])
                            # print(*Genotype_info1[0], *Genotype_info1[1:], 'Right_flank', file = out, sep = '\t')

                        elif Genotype_info1[3] in ['0|1', '1|0'] and Genotype_info1[2] in ['0|0', '1|1']:
                            NHsnp_denovo_sites.append([*Genotype_info1[0], *Genotype_info1[1:], 'Right_flank'])
                            # print(*Genotype_info1[0], *Genotype_info1[1:], 'Right_flank', file = out, sep = '\t')
                            
                            
                    snp_alleles_info_RF.pop(index_0)  

        # to keep track of chr in current itteration.
        previous_chrom = chrom

    
    vcf_in.close()
    
    # returning list of NHsnps
    return NHsnp_denovo_sites
