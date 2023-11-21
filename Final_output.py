#!/usr/bin/env python

#
def merge_results(Allele_sharing, Read_tracing):

    # we save regions for both methods and only print intersection of Allele_sharing with Read_tracing.
    PofO_final_regions           = []
    Allele_sharing_regions       = []
    Read_tracing_regions         = []

    # saving the regions to check union of it.
    AS = [Allele_sharing_regions.append(f'{info[0]}-{info[1]}')  for info in Allele_sharing]
    RT = [Read_tracing_regions.append(f'{info[0]}-{info[1]}')    for info in Read_tracing_regions]       

    # AS common with RT are filterd.
    Final_regions = set(Allele_sharing_regions) - set(Read_tracing_regions)


    # merging both results
    AS_final = [PofO_final_regions.append(info)    for info in Allele_sharing   if f'{info[0]}-{info[1]}' in Final_regions]
    RT       = [PofO_final_regions.append(info)    for info in Read_tracing_regions]


    
    return PofO_final_regions