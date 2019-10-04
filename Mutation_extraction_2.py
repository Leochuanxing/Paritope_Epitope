import copy
import numpy as np
import math
import json
import os
import pandas as pd

from pyexcel_ods import get_data, save_data

#import math
os.chdir('/home/leo/Documents/Database/Pipeline_New/Codes')
from AAC_2 import  Get_contact
from FrameConstraint import Get_consecutive
os.chdir('/home/leo/Documents/Database/Pipeline_New/Codes')
from RBFN_coverage import Distance_matrix, Design_matrix

##########################################################

from Bio import Align
aligner = Align.PairwiseAligner()
aligner.open_gap_score = -1 # it is unlikly that two sequence of the same source could be different by a gap
aligner.extend_gap_score = -1
aligner.match = 1
aligner.mismatch = -1 
aligner.mode = 'global'
#####################################################################

'''
Criteria:
    First we pick the number of mutations positions, they should be continuous in the 1_free sense.
    
    Second, we find all the opposite consecutive antigen amino acids in the 1_free sense.
             Among all those consecutive antigen amino acids with length no more than 4, we 
             choose the one with the largest contact numbers with the given mutation positions.   
'''
##################################################################################
'''
Coordinates:
    a function to extract coordinates. This one is designed for Affinity validation, and it is talored from the
    'Coordinates' in AAC_2_2 but more efficient.
Input:
    file, combined_chain_id: They are the same as in AAC_2
    mutations_pos: a list gives the mutation position
    mutaion_chain: a string gives the name of the mutation chain
'''
def Coordinates(file, combined_chain_id, mutations_pos, mutation_chain):
    # creat an empty dictionary to contain the results
    cdn = {}
    for i in combined_chain_id[0]:
        cdn['h1'+i], cdn['h2'+i], cdn['h3'+i] = [], [], []
    for i in combined_chain_id[1]:
        cdn['l1'+i], cdn['l2'+i], cdn['l3'+i] = [], [], []
    for i in combined_chain_id[2]:
        cdn[i] = []
        
    # creat a tracker dictionary, and a counter dictionary
    tracker = {}
    counter = {}
    ids = combined_chain_id[0] + combined_chain_id[1] + combined_chain_id[2]
    for i in ids:
        tracker[i] = ''
        counter[i] = -1
        
    # creat a dictionary to indicate the types of chains
    chain_type = {}
    for i in combined_chain_id[0]:
        chain_type[i] = 'H'
    for i in combined_chain_id[1]:
        chain_type[i] = 'L'
    for i in combined_chain_id[2]:
        chain_type[i] = 'A'
        
#    l_range = list(range(1000))
    l_range = []
    l_range.extend(list(range(23, 41)))
    l_range.extend(list(range(49, 64)))
    l_range.extend(list(range(89, 111)))
#        l_range = [[23, 40], [49, 63], [89, 110]]
#        h_range = [[25, 37], [50, 71], [99, 129]]
#    h_range = list(range(1000))
    h_range = []
    h_range.extend(list(range(25, 38)))
    h_range.extend(list(range(50, 72)))
    h_range.extend(list(range(99,130)))
    '''Set the l_range, h_range, and a_range'''
    if mutation_chain in combined_chain_id[0]:
        for pos in mutations_pos:
            if pos not in h_range:
                return []
        h_range = mutations_pos
        l_range = []
        a_range = list(range(10000))
    elif mutation_chain in combined_chain_id[1]:
        for pos in mutations_pos:
            if pos not in l_range:
                return []
        h_range = []
        l_range = mutations_pos
        a_range = list(range(10000))
    elif mutation_chain in combined_chain_id[2]:
        a_range = mutations_pos

    
    # extract the coordinates
    for line in file:
        if line[:6] == 'ATOM  ' and line[21] in ids:
            # update the parameters
            if tracker[line[21]]!= line[22:27]:
                counter[line[21]] += 1
                tracker[line[21]] = line[22:27]
            # extract all the parameters corresponding to line[21]
            c_type = chain_type[line[21]]
            count = counter[line[21]]
            # collect the coordinates according to c_type
            if c_type == 'H':
                if count in h_range:
                    cdn['h1'+line[21]].append([float(line[30:38]), float(line[38:46]), float(line[46:54]),
                                       count, line[17:20]])
    
            if c_type == 'L':
                if count in l_range:
                    cdn['l1'+line[21]].append([float(line[30:38]), float(line[38:46]), float(line[46:54]),
                                       count, line[17:20]])

            if c_type == 'A':
                if count in a_range:
                    cdn[line[21]].append([float(line[30:38]), float(line[38:46]), float(line[46:54]),
                                           count, line[17:20]])            
    
    return cdn
#################################################################################
      
'''
Max_contact:
    returns the nple with the most contact
Input:
    consecutive:
        a list of list of consecutives with the same length. It shouldn't be empty
    selected_contact:
        a sub set of chain_contact, selected by the position
    Ab_Ag:
        A string, takes the value of either 'Ab' or 'Ag'
Output:
    max_nple:
        a list of consecutives with the maximum contact number
'''

def Max_contact(consecutive, selected_contact, Ab_Ag):
    max_contact = -1
    max_nple = []
    if Ab_Ag == 'Ab':
        chain_pos = 2
    elif Ab_Ag == 'Ag':
        chain_pos = 1
    for nple in consecutive:
        n_contact = 0
        for fcdn in selected_contact:
            if fcdn[chain_pos] in nple:                
                n_contact += fcdn[3]
        if n_contact >= max_contact:
            max_contact = n_contact
            max_nple = nple
    return max_nple
############################################################################
'''
Order_Ag_sequence:
    Change the order of the paires, so that they match
Input:
    Ab_sequence:
        a list of numbers, give the positions of the Ab. 
        #### We always suppose the Ab_sequence is ordered from small to large and fixed###
    Ag_sequence:
        a list of numbers, gives the positions of the Ag
        #### Ag_sequence is also arranged from small to large, but it may be reversed 
             according to the criteria.#####
    contact:
        the four coordinates contact information about the amino acids between Ab_sequence and
        Ag_sequence
Output:
   Ag_sequence:  
       a list of numbers with either the same order as Ag_sequence or the reversed order.
'''
def Order_Ab_sequence(Ab_sequence, Ag_sequence, contact):
    start_para = Ab_sequence[0]
    end_para = Ab_sequence[-1]
    start_epi = Ag_sequence[0]
    end_epi = Ag_sequence[-1]
    # Use the same method as what we do in the FrameConstraint
    keep = 0
    rever = 0
    for fcdn in contact:
        if fcdn[1]==start_para and fcdn[2] == start_epi:
            keep += fcdn[3]
        elif fcdn[1] == end_para and fcdn[2] == end_epi:
            keep += fcdn[3]
        elif fcdn[1] == start_para and fcdn[2] == end_epi:
            rever += fcdn[3]
        elif fcdn[1] == end_para and fcdn[2] == start_epi:
            rever += fcdn[3]
    if rever > keep:
        Ab_sequence.sort(reverse=True)

####################################################################
'''
Select_contact_opposite:
    a function to select the contact from among all the contact
    and the all the possible opposite
Input:
    mutation_match_parameter:
        a dictionary, contains
        mutation_match_parameter['pdbid']
        mutation_match_parameter['mutation_chain']: 
            a string, gives one of the names of the mutation chains
        mutation_match_parameter['mutations_pos']:
            a list of integers, gives the positions of the mutations on the mutations_chain, it should
            in an increasing order.
        mutation_match_parameter['matched_ids']:
        mutation_match_parameter['combined_ids']:###Can we creat the combined ids according to the matched ids?
            the combined ids corresponding to the pdbid
        mutation_match_parameter['opposite_chain']:
            a string, give one of the names of the opposite chains
        mutation_match_parameter['Ab_Ag']: gives the information of whether the 
                                           mutation chain is an Ab chain or an Ag chain.
       
        cutoff: 
            a float, give the relaxed cutoff value
'''
def Select_contact_opposite(mutation_match_parameter, sequence, cutoff,\
                            moving, moving_step, moving_start, moving_end):
    # take out the parameter
    pdbid = mutation_match_parameter['pdbid']
    chain = mutation_match_parameter['mutation_chain']
    mutations_pos = mutation_match_parameter['mutations_pos']
    matched_ids = mutation_match_parameter['matched_ids']
    combined_ids = mutation_match_parameter['combined_ids']
    opposite_chain = mutation_match_parameter['opposite_chain']
    Ab_Ag = mutation_match_parameter['Ab_Ag']
    # Extract the required data
    os.chdir('/home/leo/Documents/Database/Data_Code_Publish/Structures/imgt')
    with open(pdbid+'.pdb', 'r') as file:
        cdn = Coordinates(file, combined_ids, mutations_pos, mutation_chain=chain)
    # Check if the cdn is []
    if cdn == []:
        return [],[],[]
#        pdbid_sequence = Chain_seq(file, combined_ids)
    pdbid_sequence = sequence[pdbid]
    # store the paires in paires
    possible_opposite = []
    if Ab_Ag == 'Ab':
        chain_pos = 2; opposite_chain_pos = 3; aa_pos = 1; opposite_aa_pos = 2
    elif Ab_Ag == 'Ag':
        chain_pos = 3; opposite_chain_pos = 2; aa_pos = 2; opposite_aa_pos = 1
    
    selected_contact = []; possible_opposite = []
    if not moving:    
        contact = Get_contact(cdn, matched_ids, cutoff)
        # Carry out the above process:
        # take out all the contact containing the chain_name
        '''*******This part can be modified to add the moving cut mode*******'''        
        if contact != []:
            for i in contact:
                if chain == i[0][chain_pos] and i[aa_pos] in mutations_pos \
                and i[0][opposite_chain_pos] == opposite_chain:                        
                    selected_contact.append(i)
                    possible_opposite.append(i[opposite_aa_pos]) 
                    
    if moving:
        moving_cutoff = moving_start
        while possible_opposite == [] and moving_cutoff <= moving_end:            
            contact = Get_contact(cdn, matched_ids, moving_cutoff)       
            if contact != []:
                for i in contact:
                    if chain == i[0][chain_pos] and i[aa_pos] in mutations_pos \
                    and i[0][opposite_chain_pos] == opposite_chain: 
                           
                        selected_contact.append(i)
                        possible_opposite.append(i[opposite_aa_pos]) 
            moving_cutoff += moving_step
                                    
    return selected_contact, possible_opposite, pdbid_sequence


###################################################################
'''
Paire_select:
    a function to carry out the above criteria
Input:
    mutation_match_parameter, the same as the input of Select_contact_opposite
        
Output:
    mutation_match_parameter with one more key value
    mutation_match_parameter['pairs']:
        gives the paires of the paratope and epitope of the complex with pdb id pdbid.
        
        #### The paratope are/is corresponding to the given positions of the mutations
        #### The epitope are the ones with the longest consecutive sequences in the 1_free 
        #### sense, but not more than four amino acids
'''
def Paire_select(mutation_match_parameter, sequence, mode, cutoff,\
                 moving, moving_step, moving_start, moving_end, one_to_one):

    # Extract the required information from the pdb file
    selected_contact, possible_opposite, pdbid_sequence =\
        Select_contact_opposite(mutation_match_parameter, sequence, cutoff,\
                                moving, moving_step, moving_start, moving_end)
    # Take out the values from the parameters
    chain = mutation_match_parameter['mutation_chain']
    mutations_pos = mutation_match_parameter['mutations_pos']
    opposite_chain = mutation_match_parameter['opposite_chain']
    Ab_Ag = mutation_match_parameter['Ab_Ag']
    # make sure we have got something and Get the longest possible consecutives 
    # The basic assumuption is that the more information we know about the surrounding 
    # environment, the more accurate our prediction. This is the reason why we choose
    # the longest possible consecutives.
    paires = []
    if possible_opposite != []:
        # Find all the longest possible consecutives
        possible_opposite.sort()
        '''***********************'''
        if one_to_one:
            length_list = [1]
        if not one_to_one:
            length_list = [3, 2, 1]
        '''************************'''    
        for length in length_list:
            longest_possible_consecutives = Get_consecutive(possible_opposite, length, free_type=0)
            if longest_possible_consecutives != []:
                break
        if mode == 'all':   
            for choosen_opposite_pos in longest_possible_consecutives:
                temp_paires = Sub_paire_select(mutations_pos, choosen_opposite_pos, selected_contact,\
                         Ab_Ag, pdbid_sequence, chain, opposite_chain)
                paires.extend(copy.deepcopy(temp_paires))
        elif mode == 'single':
            maxt_n = Max_contact(longest_possible_consecutives, selected_contact, Ab_Ag)
            paires = Sub_paire_select(mutations_pos, maxt_n, selected_contact,\
                         Ab_Ag, pdbid_sequence, chain, opposite_chain)
    
    # Load the results
    mutation_match_parameter['paires'] = paires
            
#####################################################################
def Sub_paire_select(mutations_pos, choosen_opposite_pos, selected_contact,\
                     Ab_Ag, pdbid_sequence, chain, opposite_chain):
    # Creat the empty container
    paires = []
    # If the length mutaions position is out of range, then the prediction is out of domain
    if len(mutations_pos) > 4:
        return paires
    
    choosen_opposite_pos.sort()
    # define a small function to change the order of the paires
    #        if len(mutations) >= 2 and len(choosen_opposite_pos)>=2:
    if len(choosen_opposite_pos)>=2:
        if Ab_Ag == 'Ab':
             Order_Ab_sequence(mutations_pos, choosen_opposite_pos, selected_contact)
        elif  Ab_Ag == 'Ag':
            Order_Ab_sequence(choosen_opposite_pos, mutations_pos, selected_contact)            
    
    # Load the amino acids to the paires according to the choosen_epitope_pos          
    original_aa = []
    for i in mutations_pos:
        original_aa.append(pdbid_sequence[chain][i])
    opposite_aa = []
    for i in choosen_opposite_pos:
        opposite_aa.append(pdbid_sequence[opposite_chain][i]) 
    # Make a deep copy to be safe       
    kept_opposite_aa = copy.deepcopy(opposite_aa)
    kept_original_aa = copy.deepcopy(original_aa)
    # Here we have to make sure the Ab amino acids is the first element of 
    # the paires and the Ag amino acids is the second element of the paires.
    if Ab_Ag == 'Ab':
        paires.append([kept_original_aa, kept_opposite_aa])
    elif Ab_Ag == 'Ag':
        paires.append([kept_opposite_aa, kept_original_aa]) 
        
    return paires
######################################################################

'''
Original_mutation_sets:
    to get the original matched parepi sets and the mutationed parepi sets
Input:
    mutation_match_parameter:
        with one more keys than the above
        
        ###mutation_match_parameter['mutations_aa']= [aa,aa]###
        
        a list of amino acids corresponding tho the mutations

Output:
    mutation_match_parameter with one more key value
    mutation_match_parameter['original_mutation_sets']:
        in the form of [[sets of original parepi],[sets of mutated parepi]]
'''
def Original_mutation_sets(mutation_match_parameter):
    # take out the values
    mutate_to_aa = mutation_match_parameter['mutate_to_aa']
    Ab_Ag = mutation_match_parameter['Ab_Ag']
    paires = mutation_match_parameter['paires']
    # Do the work
    if paires != []:
        mutated_paires = []
        if Ab_Ag == 'Ab':
            for parepi in paires:
                mutated_paires.append([mutate_to_aa, parepi[1]])
        elif Ab_Ag == 'Ag':
            for parepi in paires:
                mutated_paires.append([parepi[0], mutate_to_aa])
        # Load the results to original_mutation_sets.
        original_mutation_sets = [paires, mutated_paires]
        mutation_match_parameter['original_mutation_sets'] = original_mutation_sets
    else:
        mutation_match_parameter['original_mutation_sets'] = []


######################################################################

'''
Input:
    a list of the form 
   [ ['1hh6', 'C', 7, 'ASN'], ['', 'C', 9, 'LYS'], ['', 'C', 11, 'TYR'], ['', 'C', 15, 'SER'], \
   ['', 'D', 9, 'LYS'], ['', 'D', 11, 'TYR'], ['', 'D', 15, 'SER']]
Output:
    a list of the form 
    [['C', [7, 9, 11], ['ASN', 'LYS', 'TYR']],
 ['C', [15], ['SER']],
 ['D', [9, 11], ['LYS', 'TYR']],
 ['D', [15], ['SER']]]
'''
def Sub_parse_mutations(mutations_list_onepdb, free_type, single):
    if single:
        free_type = -1
    # Empty container for the final retrun value
    chain_pos_aa_list = []
    # Creat an empty dictionary
    chain_mutation_dict = {}
    for mutation in mutations_list_onepdb:
        if mutation[1] not in chain_mutation_dict:
            chain_mutation_dict[mutation[1]] = [mutation[2:4]]
        else:
            chain_mutation_dict[mutation[1]].append(mutation[2:4])
            
    # Find the consecutives for each chain        
    for key, value in chain_mutation_dict.items():
        value.sort(key = lambda x:x[0])
        
        # Find the consecutive pos and the consecutive aa
        chain_pos_aa_list.append([key, [value[0][0]], [value[0][1]]])
        if len(value) > 1:
            for i in range(1, len(value)):
                diff = value[i][0] - chain_pos_aa_list[-1][1][-1]
                if diff <= free_type+1:
                    chain_pos_aa_list[-1][1].append(value[i][0])
                    chain_pos_aa_list[-1][2].append(value[i][1])
                else:
                    chain_pos_aa_list.append([key, [value[i][0]], [value[i][1]]])
                    
    return chain_pos_aa_list
########################################################################
'''
Parse_mutations:
    This function is to pares the mutation information into a list of dictionaries
Input:
    mutation information from the .ods file
Output:
    list_mutaions
    a list of dictionaries with keys pdbid and values dictionries, 
    the value dictionary has keys 'mutations' and 'affinities'
**********************************************************************
The small_data will be parsed into the following form.
[[['1hh9',[ 'C', [7, 9], ['GLY', 'ARG']],['affinities', 1.00E-05, 1.00E-07, 'Kd']],
  [['1hh6', ['C', [7, 9], ['ASN', 'LYS']], [affinities', 1.00E-07, 1.00E-05], 'Kd']]]
'''
def Parse_mutations(mutation_data, free_type, single):
    # Separate the mutation complexes
    separated_mutation_complex = []
    starts = []
    ends = []
    for i in range(len(mutation_data)):
        if len(mutation_data[i][0]) == 4:
            starts.append(i)
        elif mutation_data[i][0] == 'affinities':
            ends.append(i)
    
    if len(starts) != len(ends):
        print(len(starts),' complexes \n', len(ends), ' affinities')
        return None
    
    for nth_complex in range(len(starts)):
        start = starts[nth_complex]
        end = ends[nth_complex]
        
        # Get the pdbid
        pdbid = mutation_data[start][0]
        # Get the mutation id
        mutation_id = mutation_data[start][-1]
        # get the list of mutaions
        mutations_list_onepdb = mutation_data[start:end]    
        # grouped mutations by free_type
        mutations_by_free_type_chain = Sub_parse_mutations(mutations_list_onepdb, free_type, single)
        # affinities
        affinities = mutation_data[end]
        
        separated_mutation_complex.append([pdbid, mutations_by_free_type_chain,affinities, mutation_id])
        
    return separated_mutation_complex



################################################################################
'''
Initiate_mutation_match_parameter:
    a function to initiate the mutation_match_parameter for later usage
Input:
    good_matched_ids:
        The same as we always use.
    good_combined_ids:
        As above
    one_mutation:
        [pdbid, mutations_by_free_type_chain,affinities, mutation_id]
Output:
    unit_list:
        a list of mutation_match_parameters, and this list is a basic prediction unit
'''
#math.exp(1*1000/(8.31*298))
def Initiate_mutation_match_parameter(matched_ids, combined_ids, one_mutation):
    # Take out the values from the one_mutations
    pdbid = one_mutation[0] 
    # Mutation_id
    mutation_id = one_mutation[3]
    # Find the affinity type  
    affinity_type = one_mutation[2][3]    
    # Calculate the fold change
    '''Here we would like to transform all the fold change to Kd '''
    if affinity_type == 'Kd':
        fold_change = one_mutation[2][2]/one_mutation[2][1]
    elif affinity_type == 'Ka':
        fold_change = one_mutation[2][1]/one_mutation[2][2]
    elif affinity_type == 'DDG':
        if one_mutation[2][2] == 0:
            fold_change = 1
        else:
            fold_change = math.exp(one_mutation[2][2]*1000/(8.31*298))
    else:
        print(affinity_type, pdbid)
    # Find the combined ids
    combined_ids = combined_ids[pdbid]
    # Find the matched ids
    matched_ids = matched_ids[pdbid]
    # Get the opposite chain
    # create an empty list
    unit_list = []
    # Load up the list
    for mutations in one_mutation[1]:
        mutation_chain = mutations[0]
        mutations_pos = mutations[1]
        mutate_to_aa = mutations[2]        
        opposite_chains = []
        for match in matched_ids:            
            for i in range(3):
                if match[i] == mutation_chain and i == 2:
                    Ab_Ag = 'Ag'
                    for j in range(2):
                        if match[j] not in opposite_chains and match[j] != '':
                            opposite_chains.append(match[j])                  
                elif match[i] == mutation_chain and i != 2:
                    Ab_Ag = 'Ab'
                    opposite_chains.append(match[2])

        for opposite in opposite_chains:
            # Creat an empty dictionary
            temp_mutation_match_parameter = {}
            # Load the Ab_Ag
            temp_mutation_match_parameter['Ab_Ag'] = Ab_Ag
            # Load the pdbid
            temp_mutation_match_parameter['pdbid'] = pdbid
            # Load the mutation id
            temp_mutation_match_parameter['mutation_id'] = mutation_id
            # Load the mutation chain
            temp_mutation_match_parameter['mutation_chain'] = mutation_chain
            # Load the mutations
            temp_mutation_match_parameter['mutations_pos'] = mutations_pos
            # Load the mutations_aa
            temp_mutation_match_parameter['mutate_to_aa'] = mutate_to_aa
            # Load the opposite chain 
            temp_mutation_match_parameter['opposite_chain'] = opposite
            # Load the fold change
            temp_mutation_match_parameter['fold_change'] = fold_change
            # Load the matched_ids
            temp_mutation_match_parameter['matched_ids'] = matched_ids
            # Load the combined_ids
            temp_mutation_match_parameter['combined_ids'] = combined_ids
            # Load the affinity type
            temp_mutation_match_parameter['affinity_type'] = affinity_type                        
            # Add the dictionary to the unit_list
            unit_list.append(copy.deepcopy(temp_mutation_match_parameter))
    
    return unit_list


############################################################


























