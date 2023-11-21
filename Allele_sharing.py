#!/usr/bin/env python

def Allele_sharing_method(denovo_region_list):

    Allele_sharing_regions = []

    for region in denovo_region_list:

        rg = f'{region[0]}-{region[1]}'
        
        # Trio alleles
        O_gts = region[5].split('/')
        F_gts = region[6].split('/')
        M_gts = region[7].split('/')
        Parent_allels = F_gts+M_gts


        # Here we ignore homo-alt and 
        if O_gts[0] not in Parent_allels and O_gts[1] not in Parent_allels : continue # Homozygous alternative skiped.
        if Parent_allels.count(O_gts[0]) >= 3 or Parent_allels.count(O_gts[1]) >= 3 : continue # cases like, O : 10|25   P : 10|10   M : 10|10 or 11|10
        
        # de novo allele len and non de novo allele len.
        denovo_allele     =  region[-1][0]
        non_denovo_allele =  O_gts.copy()
        non_denovo_allele.remove(str(denovo_allele))
        
        # Here we ignore case which is not possible with allele sharing. 20|21  21|23  19|21 -- non-denovo allele is having possibility of inheritance both parents.
        if non_denovo_allele[0] in F_gts and non_denovo_allele[0] in M_gts: continue 
        
        # Allele sharing method.
        parents_alleles        =  F_gts+M_gts
        parent_identification  =  ['M', 'M', 'P', 'P']

        
        PofO_determined = ''
        # Checking all the cases of Allele sharing for phasing DSTRs
        if parents_alleles.count(denovo_allele) > 0 :
            PofO_determined = parent_identification[parents_alleles.index(denovo_allele)]
            Allele_sharing_regions.append([*region[:-1]]+[denovo_allele]+[PofO_determined])

        
        elif parents_alleles.count(non_denovo_allele[0]) > 0 :
            PofO_determined = parent_identification[parents_alleles.index(non_denovo_allele[0])]
            Allele_sharing_regions.append([*region[:-1]]+[denovo_allele]+[PofO_determined])


    return Allele_sharing_regions