#!/usr/bin/env python

# required packages

import gzip
from tqdm import tqdm
from collections import Counter
import os
import numpy as np
import itertools
import pysam
from pysam import VariantFile
import ast


# Function to read the MonSTR vcf file and identify all possible de novo cases. cases like (O:12|14   F:16|18   M:20|22) de novo allele is confidently picked at the PofO level, after getting phased.
def identify_denovo(trio_vcf_file):


    with gzip.open(trio_vcf_file, 'rt')  as GangSTR_vcf:
       
        denovo_region_list = [] # let's store all the de novo regions in list of list, coz we pass it to NHsnp_identification module.

        for i,region in tqdm(enumerate(GangSTR_vcf)):
            
            if len(denovo_region_list) > 20: break
            # As it is vcf file read as .tsv we ignore hearder rows.
            if region.startswith('#'): continue
            region = region.strip().split('\t')
    
            # continuing the same as reference and failed regions.
            if region[4] == "."  : continue # if Trio family has alleles same as reference.
            if not region[-3][-4:] == region[-2][-4:] == region[-1][-4:] == 'PASS' : continue  # O, F, M info field has filter info mentioned at the end, which must be equal to "PASS"
    
            # identifying the Genotype lengths from the alternative alleles & considering ref seq as repeat lenght.
            info_field = region[7].split(';')
            end = info_field[0].split('=')[-1]
            motif = info_field[1].split('=')[-1]
            start = region[1]
            chrom = region[0]
            rep_len =  len(region[3]) #int(end) - int(start)  # in future if we want to use repeat length as "end - start" 
            # ref_len =  ['A'* rep_len] # in future if we want to use repeat length as "end - start"
            family_Alleles =  [region[3]] + region[4].split(',')
    
            # below 3 variables holding family genotypes.
            F_GTs = [int(region[-3][0]), int(region[-3][2])] # region[-3] has '0/1' or '1/0' etc, where [ [-3][0], [-3][2] ] = [0, 1] or [1, 0] 
            M_GTs = [int(region[-2][0]), int(region[-2][2])]
            O_GTs = [int(region[-1][0]), int(region[-1][2])]
    
            # getting the allele lengths for O,F,M from the alt allele reported. As alt seq is reported in alt col, based on reported genotype (0/1, 1/0, 2/1,3/0) we can directly index the length of the seq as allele length.
            F_alleles = [len(family_Alleles[F_GTs[0]]), len(family_Alleles[F_GTs[1]])]
            M_alleles = [len(family_Alleles[M_GTs[0]]), len(family_Alleles[M_GTs[1]])]
            O_alleles = [len(family_Alleles[O_GTs[0]]), len(family_Alleles[O_GTs[1]])]           
    
            
            # checking heterozygous alleles, --- tough to explain in comments ask me how this logic works. me = kiran.
            Het_Alt_Homo_Alt = set(O_alleles)-(set(F_alleles+M_alleles))
            offspring_alleles = list(set(O_alleles))[0]
            
            current_line = []
            # Checking if off-spring having any de novo allele.
            if O_alleles == M_alleles  and  len(set(O_alleles).intersection(set(F_alleles))) == 0   or  O_alleles == F_alleles  and  len(set(O_alleles).intersection(set(M_alleles))) == 0 :
                                                                
                denovo_region_list.append([chrom, start, end, rep_len, motif, f'{O_alleles[0]}/{O_alleles[1]}', f'{F_alleles[0]}/{F_alleles[1]}', f'{M_alleles[0]}/{M_alleles[1]}', list(set(O_alleles))])    

            
            elif (len(Het_Alt_Homo_Alt) > 0 and len(set(O_alleles)) == 2) or (len(set(O_alleles)) == 1 and  (F_alleles+M_alleles).count(offspring_alleles) == 1) or (len(set(O_alleles)) == 1 and (F_alleles+M_alleles).count(offspring_alleles) == 0):                
                if len(Het_Alt_Homo_Alt) == 0:
                    denovo_region_list.append([chrom, start, end, rep_len, motif, f'{O_alleles[0]}/{O_alleles[1]}', f'{F_alleles[0]}/{F_alleles[1]}', f'{M_alleles[0]}/{M_alleles[1]}', list(set(O_alleles))])
                    
                else: 
                    denovo_region_list.append([chrom, start, end, rep_len, motif, f'{O_alleles[0]}/{O_alleles[1]}', f'{F_alleles[0]}/{F_alleles[1]}', f'{M_alleles[0]}/{M_alleles[1]}', list(Het_Alt_Homo_Alt)])


        return denovo_region_list
        
