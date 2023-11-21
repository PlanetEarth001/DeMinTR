#!/usr/bin/env python

import gzip
from tqdm import tqdm
from collections import Counter
from pysam import VariantFile
import pysam
import ast



# Here are the couple of offspring-cases where Read-tracing method becomes ambigious,  (O, F, M) - ("12/12"  12/12  14/16) , ("12/14"  12/14  16/18).

def Ambiguity_identification(Determine_PofO_results):

    Read_tracing_final_regions = []

    for line in Determine_PofO_results:

        region = f'{line[0]}-{line[1]}'
        
        # info we needed.
        O_alleles = line[5]
        F_alleles = line[6].split('/')
        M_alleles = line[7].split('/')
        denovo_allele = line[8][0]
        snp_info = line[10].split(';')
        parents_gts = {'P':F_alleles, 'M':M_alleles}
        snp_pofo = {snp_info[0].split(':')[0] : snp_info[0].split(':')[1], snp_info[1].split(':')[0] : snp_info[1].split(':')[1]}
        determined_pofo = snp_pofo[line[-4]]
        determined_pofo_snp_allele = line[-4]

        
        # Offspring alleles being same. 12|12...
        if O_alleles.split('/')[0] == O_alleles.split('/')[1]:

            # Allele sharing method to phase.
            parents_alleles = F_alleles+M_alleles
            parent_identification = ['P', 'P', 'M', 'M']
            
            PofO_determined = ''
            # Checking all the cases of Allele sharing for phasing DSTRs
            if parents_alleles.count(str(denovo_allele)) > 0 :
                PofO_determined = parent_identification[parents_alleles.index(str(denovo_allele))]

                # Checking if snp-inherited parent is same as allele-sharing-inherited parent.
                if determined_pofo == PofO_determined :         
                    snp_allele = list(snp_pofo.keys())
                    snp_allele.remove(determined_pofo_snp_allele)
                    Read_tracing_final_regions.append([*line]+[snp_pofo[snp_allele[0]]])
                elif determined_pofo != PofO_determined :
                    Read_tracing_final_regions.append([*line]+[snp_pofo[determined_pofo_snp_allele]])
            else: Read_tracing_final_regions.append([*line]+[het_alt])


        
        # Offspring alleles being diff, 12|14  12|14  12|13
        elif O_alleles in [line[7], line[6]]:
            associated_len = line[-5]

            # checking if the associated allele len is in which parent.
            if associated_len in parents_gts[determined_pofo]:
                denovo_allele_case2 = eval(line[8])
                denovo_allele_case2.remove(int(associated_len))
                gts = list(parents_gts.keys())
                gts.remove(determined_pofo)
                Read_tracing_final_regions.append([*line]+[gts[0]])
            else: 
                Read_tracing_final_regions.append([*line]+[determined_pofo])

        # regions which doesnt fall under  above 2 cases.
        else : Read_tracing_final_regions.append([*line])


    return Read_tracing_final_regions