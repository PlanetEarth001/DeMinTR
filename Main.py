#!/usr/bin/env python
# 

# importing packages required for main function.
import argparse

# importing required files as modules.
from denovo_utils        import identify_denovo
from NHsnps_utils        import Proximal_snps
from Determination_PofO  import Determine_PofO 
from Ambiguity_utils     import Ambiguity_identification
from Allele_sharing      import Allele_sharing_method
from Final_output        import merge_results


# Point to notice, this file is main file and "parse_args" function will run always & __name__ == "__main__" considers as main function.
def parse_args():
    parser = argparse.ArgumentParser( prog='Eta', description='Eta is a tool to identify STR de novo mutations and determine their PofO', epilog='')

    # adding arguments to instance of class, which displays at command line.    
    parser.add_argument('--Joint_Trio_vcf',    type=str,  dest='Joint_Trio_vcf',         required=True,  help= 'Input : Joint Trio STR VCF file')
    parser.add_argument('--Joint_vcf_loc',     type=str,  dest='Joint_vcf_loc',          required=True,  help= 'Input : Joint SNP|INDEL vcf file')
    parser.add_argument('--Offspring_file_id', type=str,  dest='Offspring_file_id',      required=True,  help= 'Input : Offspring file id')
    parser.add_argument('--Paternal_file_id',  type=str,  dest='Paternal_file_id',       required=True,  help= 'Input : Paternal file id')
    parser.add_argument('--Maternal_file_id',  type=str,  dest='Maternal_file_id',       required=True,  help= 'Input : Maternal file id')
    parser.add_argument('--Alignment_file',    type=str,  dest='Alignment_file',         required=True,  help= 'Input : Maternal file id')
    # filter arguments need to updated
    # filter: fragment size, number of reads, min read quality, flank disnt, Allele drop out, additional info feild
    
    
    args = parser.parse_args()

    return args




# Main function
if __name__ == "__main__":   # python interpreter consideres this as a main function.

    args = parse_args()      # args object will be usefull to retrive parameters given at user-end from command line. 


    # STEP:1 -------Passing GangSTR vcf file to denovo_utils and storing list of denovo regions in denovo_regions_list - var. (Identifying de novo)
    trio_input_vcf        =   args.Joint_Trio_vcf    
    denovo_regions_list   =   identify_denovo(trio_input_vcf)


    # STEP:2 -------Here we pass list of denovo regions to NHsnps finding module.  (Identifying NHsnps for each possible de novo site)
    joint_vcf_location    =   args.Joint_vcf_loc
    O_file_id             =   args.Offspring_file_id
    F_file_id             =   args.Paternal_file_id
    M_file_id             =   args.Maternal_file_id
    regions_with_NHsnps   =   Proximal_snps(joint_vcf_location, O_file_id, F_file_id, M_file_id, denovo_regions_list)


    # STEP:3 --------Passing NHsnps sites to PofO module.
    Alignment_file        =   args.Alignment_file
    Determine_PofO_results=   Determine_PofO(regions_with_NHsnps, Alignment_file)


    # STEP:4 --------By step-3 we generated PofO result for all possible de novo sites. "but" 
                     # cases like 12|12 or 12|14 offspring GTs creates ambiguity while confirming the de novo allele. Step-4 identifies the ambiguity cases and clears them.
    Read_tracing_result   =   Ambiguity_identification(Determine_PofO_results)


    # STEP:5 --------In this step we apply Allele sharing method.
    Allele_sharing_results=   Allele_sharing_method(denovo_regions_list)


    # STEP:6 --------Finally we got 2 methods regions (Allele sharing & Read tracing), this module intersects "Read tracing" with "Allele sharing" and prints the output.
    Final_output_print    =   merge_results(Allele_sharing_results, Read_tracing_result)    

    # let's go print the output.
    for i in Final_output_print:
        print(i)

