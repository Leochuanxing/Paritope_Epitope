import copy
import numpy as np
import math
import json
import os
import xlrd
import pandas as pd
import random
from pyexcel_ods import get_data, save_data
from scipy.stats.stats import pearsonr
from pydoc import help
#import math
os.chdir('/home/leo/Documents/Database/Data_Code_Publish/Codes')
from AAC_2 import  Get_contact
from FrameConstraint import Get_consecutive
from RBFN_coverage import Distance_matrix, Design_matrix
from matplotlib import pyplot as plt

##########################################################

from Bio import Align
aligner = Align.PairwiseAligner()
##################################
'''Here the open_gap score and the extened gap score is consistant with the best scores
    selected in the cross validation'''
aligner.open_gap_score = -1 
aligner.extend_gap_score = -1
#######################################
aligner.match = 1
aligner.mismatch = -1 
aligner.mode = 'global'
#####################################################################
#
#'''
#Criteria:
#    First we pick the number of mutations positions, they should be continuous in the 1_free sense.
#    
#    Second, we find all the opposite consecutive antigen amino acids in the 1_free sense.
#             Among all those consecutive antigen amino acids with length no more than 4, we 
#             choose the one with the largest contact numbers with the given mutation positions.   
#'''
##################################################################################
'''
Coordinates:
    a function to extract coordinates. This one is designed for Affinity validation, and it is talored from the
    'Coordinates' in AAC_2_2 but more efficient.
Input:
    file, combined_chain_id: They are the same as in AAC_2
    Ab_Ag: the mutations chain
    mutations_pos: a list gives the mutation position
    mutaion_chain: a string gives the name of the mutation chaing
'''
def Coordinates(file, combined_chain_id, Ab_Ag, mutations_pos, mutation_chain):
    # creat an empty dictionary to contain the results
    cdn = {}
    for i in combined_chain_id[0]:
        cdn['h1'+i], cdn['h2'+i], cdn['h3'+i] = [], [], []
    for i in combined_chain_id[1]:
        cdn['l1'+i], cdn['l2'+i], cdn['l3'+i] = [], [], []
    for i in combined_chain_id[2]:
        cdn[i] = []
        
    # creat a tracker dictionary, and a counter dictionary
    tracker = {}; counter = {}
    ids = combined_chain_id[0] + combined_chain_id[1] + combined_chain_id[2]
    for i in ids:
        tracker[i] = ''; counter[i] = -1
        
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
        a sub set of chian_contact, selected by the position
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
        chain_pos = 1
    elif Ab_Ag == 'Ag':
        chain_pos = 2
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
    with open(pdbid+'.pdb', 'r') as file:
        cdn = Coordinates(file, combined_ids, Ab_Ag, mutations_pos, mutation_chain=chain)
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
    
    selected_contact = []; possible_opposite = []; equal_mutations = [] 
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
                    equal_mutations.append(i[aa_pos])
                    
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
                        equal_mutations.append(i[aa_pos])
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
 
#######################################################################
'''
Prediction_RBFN_coverage:
    a function to calculate the predicted values of the testing_set
    
Input:
    testing_set:
        a list of the given Ab_aa and Ag_aa paires in the form of 
       [[[ASP, LYS], [ALA, ARG, GLU]],......]
       
    wd_results:
        gives the working directory of the results

Output:
    predictions:
        an arrey in the shape of (len(testing_data), a), gives the predicted values
        of the paires.
'''
def Predition_RBFN_coverage(train_results, testing_set):
#testing_set = testing5
    key = str(len(testing_set[0][0]))+'_'+str(len(testing_set[0][1]))+'_0_0_1_2_1perchain' 
    '''********************Pay attention to this scale factor****************'''
#    scale_factor = len(testing_set[0][0]) * len(testing_set[0][1])
    
    coeff = train_results[key]['coefficients']
    coeff = np.reshape(np.array(coeff), (-1,1))
    centers_selected = train_results[key]['centers']
    # Calculate the distance matrix
    distance_matrix = Distance_matrix(testing_set, centers_selected, square = False)

    linear_coeff = np.reshape(coeff, (-1,1))

    radius_coeff = np.ones((distance_matrix.shape[1], 1))
    testing_design_matrix = Design_matrix(distance_matrix, radius_coeff, basis_function = 'Gaussian')
    testing_design_matrix = np.hstack((testing_design_matrix, np.ones((len(testing_set), 1))))
    
    predictions = testing_design_matrix.dot(linear_coeff)
    '''******************************'''
#    predictions *= scale_factor
    
    return predictions
  
##################################################################
'''
Compare:
    a function to compare the predicted values of the original data set and
    the mutated data set
Input:
     
     wd_results: The same as above
     wd_negative_samples:The same as above
     mutation_match_parameter['fold_change']: a float number, give the fold
     of the Kd changes
     model: a string, it takes the value of either 'RBFN' or 'Logistic'
     if it takes the value of 'RBFN', we use the RBFN model to make prediction
     if it takes the value of 'Logistic', we use the Logistic regression model to make prediction.
Output:
    compare:
        a list of float numbers, [original_pred_sum, mutation_pred_sum]
'''

def Compare(mutation_match_parameter, test_coverage_RBFN_results):
    # Take out the values
#    fold_change = mutation_match_parameter['fold_change']
    original_mutation_sets = mutation_match_parameter['original_mutation_sets']
    # Get the original sets and the mutated sets
    if original_mutation_sets != []:
        original_sets = original_mutation_sets[0]
        mutated_sets = original_mutation_sets[1]
        # Make predictions for each set
        original_pred = Predition_RBFN_coverage(test_coverage_RBFN_results,original_sets)
        mutated_pred = Predition_RBFN_coverage(test_coverage_RBFN_results, mutated_sets)
        # Put together
        '''************Should we use sum or average********************'''
        original_mutation_score = [np.average(original_pred), np.average(mutated_pred)]
    else:
        return 'Empty Match'        
        
    return original_mutation_score

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
    This function is to pares the mutation information into a list of dictionaies
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

#########################################################################              
'''
Mutation_list_prediction:
    to make prediction for a mutation list
input:
    one_list: 
        a list containing a set of mutations corresponding to measured affinity
Output:
    scores:
        either[original_score, mutation_score]
        or ['Empty Match']
'''

def Mutation_list_prediction(one_list, sequence,  test_coverage_RBFN_results, mode, cutoff,\
                             moving, moving_step, moving_start, moving_end, one_to_one):  

    # Select the the paires and Generate the original and mutation sets
#    list_parameters = []
    for mutation_match_parameter in one_list:
        Paire_select(mutation_match_parameter,sequence, mode, cutoff,\
                     moving, moving_step, moving_start, moving_end, one_to_one, )
        Original_mutation_sets(mutation_match_parameter)
    # Make prediction
    # When we make prediction, we simply assume each paratope-epitope paire contribute equally to 
    # the final affinity
    all_original_mutation_score= []
    for parameter in one_list:
        original_mutation_score = Compare( parameter, test_coverage_RBFN_results)
        # Load the results to scores
        all_original_mutation_score.append(original_mutation_score)

    # Check if there are empty match, if so this prediction doesn't count 
    for j in all_original_mutation_score:
        if j == 'Empty Match':
            # If there is one empty match, the whole prediction ends and return the value 'Empty Match'
            return [j]

    # Add up the scores
    scores_original=0
    scores_mutation = 0
    length = len(all_original_mutation_score)
    for original_mutaion in all_original_mutation_score:
        '''*********Try not to devide the length***************'''
#        scores_original += original_mutaion[0]
#        scores_mutation += original_mutaion[1]
        scores_original += original_mutaion[0]/length
        scores_mutation += original_mutaion[1]/length
    
    return [scores_original, scores_mutation]

######################################################################
'''
Test the above functions
'''
#data = small_data
def Predict(data, matched_ids, combined_ids, sequence, test_coverage_RBFN_results, mode, cutoff, free_type, single,\
            moving, moving_step, moving_start, moving_end, one_to_one):
    mu_list = Parse_mutations(data, free_type, single)
    results = []

    for mutation in mu_list: 
        sub_result = []
        one_list = Initiate_mutation_match_parameter(matched_ids, combined_ids, mutation)
        scores = Mutation_list_prediction(one_list, sequence,  test_coverage_RBFN_results, mode, cutoff,\
                                          moving, moving_step, moving_start, moving_end, one_to_one)
        '''**********'''
        sub_result.append(mutation)
        for para in one_list:
            affinity_type =  para['affinity_type']
            fold_change = para['fold_change']
            sub_result.append(para['original_mutation_sets'])
        if scores != ['Empty Match']:
            r_w = Right_or_wrong(fold_change, affinity_type, scores[0], scores[1])
            scores.extend([fold_change,r_w])
        
        sub_result.append(scores)
        results.append(copy.deepcopy(sub_result))
        
    return results  
#results = Predict(data) 
#results
############################################################

'''
Right_or_wrong:
    a function to tell whether our  prediction is right or wrong
Input:
    fold: a float, shows how many times the affinity changes
    affinity_type: string, either 'Ka' or 'Kd'
    prediction1: a float, give the score of the first prediction
    predictin2: a float, give the score of the second prediction
Output:
    right_or_wrong:
        a boolean, if the prediction is correct, it takes the values of True
                   if the prediction is incorrect, it takes the values of False
'''
def Right_or_wrong(fold, affinity_type, prediction1, prediction2):
    # In the case when the affinity_type is Kd
    right_or_wrong = True
#    if affinity_type == 'Kd':
    if fold < 1 and prediction1 < prediction2:
        right_or_wrong = True
    elif fold > 1 and prediction1 > prediction2:
        right_or_wrong = True
    elif fold == 1:
        right_or_wrong = 'Undecided'
    else:
        right_or_wrong = False
    
    return right_or_wrong
########################################################################

'''
Count the number of predictions and the number of correct predictions
'''
def Count_correct(results, fold_change_cut = 1):
    total_prediction = 0
    correct_prediction = 0
    
    for res in results:
        # Take out the values
        if res[-1] != ['Empty Match']:
            fold_change = res[-1][2]
            if fold_change > fold_change_cut or fold_change <1/fold_change_cut:
                if res[-1][3] != 'Undecided':
                    total_prediction += 1
                    correct_prediction += res[-1][3]
            
    return total_prediction, correct_prediction
###############################################################################
def ROC_AUC(results_all, fold_cut):
    res_all = []
    for results in results_all:
        
        if results[-1][0]!= 'Empty Match':
            res_all.append([results[-1][2], results[-1][1]-results[-1][0]])
    res_all.sort(key = lambda x:x[1], reverse = True)
    TPR = []
    FPR = []
    total_positive = 0
    total_negative = 0
    # Count the total_positive and total_negative
    for res in res_all:
        if res[0] < 1/fold_cut:
            total_positive += 1
        elif res[0] > fold_cut:
            total_negative += 1
    # Calculate the TPR and FPR
    tp = 0; fp = 0
    Dscore = []; Dfold = []
    for res in res_all:
        if res[0] < 1/fold_cut:
            tp += 1
            TPR.append(tp/total_positive)
            FPR.append(fp/total_negative)
            Dscore.append(res[1])
            Dfold.append(res[0])
        elif res[0] > fold_cut:
            fp += 1
            TPR.append(tp/total_positive)
            FPR.append(fp/total_negative)
            Dscore.append(res[1])
            Dfold.append(res[0])
    # Calculate AUC
    AUC = 0
    for i in range(1, len(TPR)):
        AUC += TPR[i] * (FPR[i]- FPR[i-1])
    
    return TPR, FPR, AUC, Dscore, Dfold
######################################################################
class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)

###########################################################################

def My_pred(single_mutation = True, binary = True, one_to_one = False):
    os.chdir('/home/leo/Documents/Database/Data_Code_Publish/Structures')
    with open('matched_ids', 'r') as f:
        matched_ids = json.load(f)
    with open('combined_ids', 'r') as f:
        combined_ids = json.load(f)
    with open('sequence', 'r') as f:
        sequence = json.load(f)
                
    os.chdir('/home/leo/Documents/Database/Data_Code_Publish/Mutations')
#    data_raw = get_data('multimutations.ods')
    data_raw = get_data('Keating.ods')
    keys = data_raw.keys()
    
    res_dict = {}
    mod='single'; mov = True  
    
    # Use different coverage model, binary or numerical
    os.chdir('/home/leo/Documents/Database/Data_Code_Publish/Codes/Results')
    with open('train_results', 'r') as f:
        train_results = json.load(f)
#            test_coverage_RBFN_results = json.load(f)
    if binary:
        test_coverage_RBFN_results = train_results['cross_binary_Gaussian_']
    else:
        test_coverage_RBFN_results = train_results['cross_numerical_Gaussian_']
        
    sufix = 'single_'+str(single_mutation) +'_'+'mode_'+mod+'_moving_'+str(mov)+'_binary_'+str(binary)

    # Store the results
    results_dict = {}
    if not mov:
        cut_range = [5.5]
    if mov:
        cut_range = [0]
    # For different cut   
    for cut in cut_range:
#        auc_cut_dict[str(cut)] = []
#        fpr_tpr_dict[str(cut)] = []
        results = []
        for key in keys:
            print('working on:  ' + key)
            data = []
            for d in data_raw[key]:
                if d != []:
                    data.append(d)
    
            os.chdir('/home/leo/Documents/Database/Pipeline_New/All with peptide 5+ resolution 4A/structure')
            
            results.extend(Predict(data, matched_ids, combined_ids, sequence,\
                              test_coverage_RBFN_results, mode = mod, cutoff=cut,\
                              free_type = 0, single=single_mutation,\
                              moving=mov, moving_step=0.25, moving_start=3,\
                              moving_end=8, one_to_one = one_to_one))
        # Load the results
        results_dict[str(cut)] = copy.deepcopy(results)

    res_dict[sufix] = copy.deepcopy(results_dict)
                    
    # Save the results              
    affinity_results = {}
    affinity_results['results_dict'] = res_dict
#    affinity_results['auc_dict'] = auc_dict
#    affinity_results['fpr_tpr_dict'] = fpr_tpr_dict
#    affinity_results['DDG_cut'] = DDG
    affinity_results
    
    return affinity_results

###############################################################################
'''
*******************Explanation of the Predict parameters****************
mode:
    a string, it takes the value of 'single' or 'all', 'single' means if there are multiple 
    opposite amino acids with equal length, the one with the largest contact is selected, while
    'all' means all of them are considered
cutoff:
    the relaxed cutoff distance, it only works when the moving takes the value of False
free_type:
    takes the value of either 0 or 1
single:
    it takes the value of either True or False. True, means the mutations are considered
    as single amino acid mutations, and combine all those mutations together. 
    False means the mutations are condidered together.
    For example, if the mutaion positions are [1,2,5], in the case when single is true,
    they are considered as mutations [1], [2], [5]. If the single takes the value False, 
    the mutations will be considered as [1,2] and [5]
moving:
    it takes the value of either True or False
    If it takes the value of True, means when finding the opposite amino acids of the mutation,
    we start with the cutoff moving_start, and gradually increase the cuoff by the moving_step, until 
    we find the opposite amino acids or the cutoff reaches the moving_end
one_to_one:
    one_to_one, takes the value of True or False, if it is true, means we will use
    the match-type (1,1) to make predictions only. In this case, the single should take 
    the value of True
'''
#with open('affinity_results_one_to_one', 'r') as f:
#    affinity_results = json.load(f)
##############################################################################

'''
THis function is to select only the predicted results in the much simpler version
The returned results contains the mutation id, DDG and predicted score change.
'''
def Cleanse_the_results(results):
    cleansed_results = []
    for key, value in results['results_dict'].items():
        result_list = value['0']
        for res in result_list:
            if res[-1][0] !='Empty Match':
                one = [res[0][0], res[0][-2][2], res[-1][1]-res[-1][0] ,res[0][-1]]
                cleansed_results.append(copy.deepcopy(one))
    return cleansed_results



###############################################################################
def AUC_under_cut(cleansed, cut_lower, cut_upper):
    tpr = []
    fpr = []
    # Select the predictions under the given cut
    selected = []
    for re in cleansed:
        if abs(re[1]) > cut_lower and abs(re[1]) < cut_upper:
            selected.append(re)
    # Calculate the number of positvie observations and the number 
    # of negative observations
    n_positive = 0; n_negative = 0
    for pr in selected:
        if pr[1] > 0:
            n_positive += 1
        elif pr[1] < 0 :
            n_negative += 1
    if n_positive*n_negative ==0:
        return 0, [], [], selected
    # Sort the prediction
    selected.sort(key = lambda x:x[2])
    # Calculate tpr and fpr
    n=0;p=0
    for pred in selected:
        if pred[1] > 0:
            p += 1
        elif pred[1] < 0:
            n += 1
        tpr.append(p/n_positive)
        fpr.append(n/n_negative)
    # Calculate the auc
    auc = 0
    for i in range(1, len(tpr)):
        auc += tpr[i]*(fpr[i]-fpr[i-1])
    
    return auc, fpr, tpr, selected

##############################################################################
'''
This function is to connect the name of the mutations with the ids of the mutations
'''
def Name_index(sheet, matched_ids):
    name_index =[]
    # Lets collect the mutation ids first
    for i in range(1,1015):
        pdbid = sheet.cell_value(i,0).lower()        
        mut_set = (sheet.cell_value(i,4))
        if pdbid in matched_ids:
            name_index.append([pdbid, mut_set, i])
    return name_index
############################################################################
'''
This function returns the same mutations as the selected, but predicted by other method
'''
def Set_prediction_set_the_same(selected, scores_bASA, name_index):
    row_ind = []
    n_row, n_col = scores_bASA.shape
    for i in range(n_row):
        pdbid = scores_bASA['#pdb'][i]
        mutation = scores_bASA['mutation'][i]
        ind = 'NA'
        for n_i in name_index:        
            if n_i[0].upper() == pdbid and n_i[1].upper() == mutation:
                ind = n_i[2]
                break
        row_ind.append(ind)
    bASA_selected = []
    for sel in selected:
        for j in range(len(row_ind)):
            if row_ind[j] == sel[3]:
                bASA_selected.append([scores_bASA['#pdb'][j],\
                                                 scores_bASA['expt'][j], scores_bASA['comp'][j]])
                # The data provided by Keating's PostDoc sometimes has two identical mutations but
                # with different computed scores. If so, we take the first one.
                break
    return bASA_selected
################################################################################
def Sub_Auc_other(other_selected):
    tpr = [];fpr=[]
    n_positive = 0; n_negative = 0
    for pr in other_selected:
        if pr[1] > 0:
            n_positive += 1
        elif pr[1] < 0 :
            n_negative += 1
    if n_positive*n_negative == 0:
        return 0, [], []
    # Sort the prediction
    other_selected.sort(key = lambda x:x[2], reverse = True)
    # Calculate tpr and fpr
    n=0;p=0
    for pred in other_selected:
        if pred[1] > 0:
            p += 1
        elif pred[1] < 0:
            n += 1
        tpr.append(p/n_positive)
        fpr.append(n/n_negative)
    # Calculate the auc
    auc = 0
    for i in range(1, len(tpr)):
        auc += tpr[i]*(fpr[i]-fpr[i-1])
    
    return auc, fpr, tpr
#############################################################
'''
This function is to compare keating's results with my results
'''
def AUC_other(selected):

    # Get the original mutations
    loc = '/home/leo/Documents/Database/Pipeline_New/Mutations/Mutation.xlsx'
    wb = xlrd.open_workbook(loc) 
    sheet = wb.sheet_by_index(0)
    
    os.chdir('/home/leo/Documents/Database/Data_Code_Publish/Structures')
    with open('matched_ids','r') as f:
        matched_ids = json.load(f)
    name_index = Name_index(sheet, matched_ids)
        
    os.chdir('/home/leo/Documents/Database/Pipeline_New/Mutations/ABbind_compt_data/ABbind_compt_data')
    scores_bASA = pd.read_csv('bASA_all_side.dat', sep = ';')
    scores_dDFIRE = pd.read_csv('dDfire.dat', sep = ';')
    scores_DFIRE = pd.read_csv('dfire.dat', sep = ';')
    scores_Design_studio = pd.read_csv('d_studio_expt_comp.dat', sep = ';')
    scores_Fold_X = pd.read_csv('foldX_IE.dat', sep = ';')
    scores_Rosetta = pd.read_csv('rosetta_expt_comp.dat', sep = ';')
    scores_Statium = pd.read_csv('statium_expt_comp.dat', sep = ';')
    
    scores_list = [scores_bASA, scores_dDFIRE, scores_DFIRE, scores_Design_studio,\
                   scores_Fold_X, scores_Rosetta, scores_Statium ]
    keys = ['bASA', 'dDFIRE', 'DFIRE', 'Design_studio', 'FoldX', 'Rosetta', 'Statium']
    other_pred_dict = {}
    for i in range(len(keys)):
        other_pred_dict[keys[i]] = {}    
        other_selected = Set_prediction_set_the_same(selected, scores_list[i], name_index)
        other_auc,other_fpr, other_tpr = Sub_Auc_other(other_selected)
        # Load up
        other_pred_dict[keys[i]]['other_selected'] = copy.deepcopy(other_selected)
        other_pred_dict[keys[i]]['other_auc'] = copy.deepcopy(other_auc)
        other_pred_dict[keys[i]]['other_fpr'] = copy.deepcopy(other_fpr)
        other_pred_dict[keys[i]]['other_tpr'] = copy.deepcopy(other_tpr)
        
    return other_pred_dict

#################################################################################
def Show_all_auc(my_auc, other_pred_dict):
    auc_all = []
    auc_all.append(['My', my_auc])
    keys = ['bASA', 'dDFIRE', 'DFIRE', 'Design_studio', 'FoldX', 'Rosetta', 'Statium']
    for key in keys:
        auc_all.append([key, other_pred_dict[key]['other_auc']])
    return auc_all
###############################################################################
def Analyze_results(affinity_results, cut_lower, cut_upper):
    cleansed_results = Cleanse_the_results(affinity_results)           
    auc, fpr, tpr, selected = AUC_under_cut(cleansed_results, cut_lower, cut_upper)
    other_pred_dict = AUC_other(selected)
    auc_all = Show_all_auc(auc, other_pred_dict)
    # my pred dict
    my_pred_dict = {}
    my_pred_dict['auc'] = auc
    my_pred_dict['fpr'] = fpr
    my_pred_dict['tpr'] = tpr
    my_pred_dict['selected'] = selected
    
    return other_pred_dict, my_pred_dict, auc_all
#############################################################################
def Bootstrap_my_AUC(my_pred_dict):
    selected = my_pred_dict['selected']
    boot_auc = []
    for i in range(10000):
        bootstrap_samples = random.choices(selected, k=len(selected))
        auc, _, _, _ = AUC_under_cut(bootstrap_samples, cut_lower = 0, cut_upper = 1000)
        boot_auc.append(auc)
    boot_auc.sort()
    CI_95 = [boot_auc[250], boot_auc[9750]]
    my_pred_dict['CI_95'] = CI_95    
#    return CI_95
#############################################################################
def Bootstrap_other_AUC(other_pred_dict):
    for key, value in other_pred_dict.items():
        other_selected = value['other_selected']
        other_auc_list = []
        print('Bootstraping '+key)
        for i in range(10_000):
            bootstrap_other = random.choices(other_selected, k=len(other_selected))
            auc, _,_ = Sub_Auc_other(bootstrap_other)
            other_auc_list.append(auc)
        other_auc_list.sort()
        other_CI_95 = [other_auc_list[250], other_auc_list[9750]]
        value['CI_95'] = copy.deepcopy(other_CI_95)
            
##############################################################################
def Prediction_on_my_data():
    affinity_results = My_pred(single_mutation = True, binary = True, one_to_one=True) 
    prediction_on_my_data ={}
    results = affinity_results['results_dict']\
                    ['single_True_mode_single_moving_True_binary_True']['0']
    clean_results = [res[-1] for res in results]
    cleaner_results = []
    for c in clean_results:
        if c[0] != 'Empty Match':
            cleaner_results.append(c)
    cleaner_results
    cleanest_ressults = []
    for cl in cleaner_results:
        cleanest_ressults.append(copy.deepcopy(['',math.log(cl[2])*(298*8.31)/1000, cl[1]-cl[0]]))
    boundaries = [[0,1000], [0,0.5],[0.5, 1000], [1, 1000]]
    for bound in boundaries:
        auc, fpr, tpr, selected = AUC_under_cut(cleanest_ressults, cut_lower=bound[0], cut_upper=bound[1])
        # Bootstrap
        boot_auc = []
        for i in range(10_000):
            bootstrap_samples = random.choices(selected, k=len(selected))
            auc, _, _, _ = AUC_under_cut(bootstrap_samples, cut_lower = 0, cut_upper = 1000)
            boot_auc.append(auc)
        boot_auc.sort()
        CI_95 = [boot_auc[250], boot_auc[9750]]
        prediction_on_my_data[str(bound[0])+'_'+str(bound[1])]={}
        prediction_on_my_data[str(bound[0])+'_'+str(bound[1])]['auc'] = auc
        prediction_on_my_data[str(bound[0])+'_'+str(bound[1])]['CI_95'] = CI_95
        prediction_on_my_data[str(bound[0])+'_'+str(bound[1])]['fpr'] = fpr
        prediction_on_my_data[str(bound[0])+'_'+str(bound[1])]['tpr'] = tpr
        
    return prediction_on_my_data
#prediction_on_my_data_numerical = Prediction_on_my_data()
#prediction_on_my_data_numerical['0_1000']['auc']
#for key, value in prediction_on_my_data.items():
#    print(value['auc'], value['CI_95'], len(value['tpr']))
#os.chdir('/home/leo/Documents/Database/Pipeline_New/Codes/Results')
#with open('prediction_on_my_data_numerical', 'w') as f:
#    json.dump(prediction_on_my_data_numerical, f)
##############################################################################
 

         
if __name__=='__main__': 
    all_pred_dict = {}
    affinity_results = My_pred(single_mutation = True, binary = True, one_to_one=True)
    all_pred_dict['affinity_results'] = affinity_results
    boundaries = [[0,1000], [0,0.5],[0.5, 1000], [1, 1000]]
    for bound in boundaries:                 
        other_pred_dict, my_pred_dict, auc_all = \
                                Analyze_results(affinity_results,\
                                                cut_lower = bound[0], cut_upper = bound[1])
        print('Boorstraping my AUC')
        Bootstrap_my_AUC(my_pred_dict)

        Bootstrap_other_AUC(other_pred_dict)
        all_pred_dict[str(bound[0])+'_'+str(bound[1])] = {}
        all_pred_dict[str(bound[0])+'_'+str(bound[1])]['other_pred_dict'] = other_pred_dict
        all_pred_dict[str(bound[0])+'_'+str(bound[1])]['my_pred_dict'] = my_pred_dict
        all_pred_dict[str(bound[0])+'_'+str(bound[1])]['auc_all'] = auc_all

# Save the results
#os.chdir('/home/leo/Documents/Database/Pipeline_New/Complexes/Results')
#with open('one_to_one_affinity', 'w') as f:
#    json.dump(all_pred_dict, f, cls=NumpyEncoder)
#def Clean_spm():
#    os.chdir('/home/leo/Documents/Database/Pipeline_New/Results')
#    with open('one_to_one_affinity', 'r') as f:
#        one_to_one_affinity = json.load(f)
#    
#    one_to_one_affinity.keys()
#    one_to_one_affinity['affinity_results']
#    
#    clean_results = Cleanse_the_results(one_to_one_affinity['affinity_results'])
#    
#    # Select the spm    
#    os.chdir('/home/leo/Documents/Database/Pipeline_New/Codes/Results')
#    data_raw = get_data('Keating.ods')
#    keys = data_raw.keys()
#    data_list = []
#    for key in keys:
#        data_list.extend(data_raw[key])
#    # get the spm id
#    spm_id = []
#    
#    for i in range(len(data_list)-1):
#        if len(data_list[i][0])==4 and data_list[i+1][0] == 'affinities':
#            spm_id.append(data_list[i][5])
#    # Select from the cleased
#    clean_spm=[]
#    for cl in clean_results:
#        if cl[3] in spm_id:
#            clean_spm.append(cl)
#            
#    return clean_spm
#############################################################################

#def Concentration(clean_results, top_percent_pred):
#    clean_results_sorted = copy.deepcopy(clean_results)
#    clean_results_sorted.sort(key = lambda x:x[2], reverse = True)
#    top_pred = clean_results_sorted[:math.floor(top_percent_pred*len(clean_results_sorted))+1]
#    negative = 0
#    for t in top_pred:
#        if t[1] <0:
#            negative += 1
#    negative /= len(top_pred)
#    
#    return negative
    
#negative = Concentration(clean_results, top_percent_pred=0.05)
#negative
##################################################################################
'''Write the prediction results in a .dat file'''
def Write_pred_results_in_dat():
    os.chdir('/home/leo/Documents/Database/Pipeline_New/Results')
    with open('one_to_one_affinity', 'r') as f:
        one_to_one_affinity = json.load(f)
    
    clean_results = Cleanse_the_results(one_to_one_affinity['affinity_results'])
    
    loc = '/home/leo/Documents/Database/Pipeline_New/Mutations/Mutation.xlsx'
    wb = xlrd.open_workbook(loc) 
    sheet = wb.sheet_by_index(0)
    
    os.chdir('/home/leo/Documents/Database/Pipeline_New/All with peptide 5+ resolution 4A')
    with open('matched_ids','r') as f:
        matched_ids = json.load(f)
    name_index = Name_index(sheet, matched_ids)

    
    os.chdir('/home/leo/Documents/Database/Pipeline_New/Mutations/ABbind_compt_data/ABbind_compt_data')
    scores_bASA = pd.read_csv('CM.dat', sep = ',')
        #row_ind = []
    n_row, n_col = scores_bASA.shape
    for i in range(n_row):
        pdbid = scores_bASA['#pdb'][i]
        mutation = scores_bASA['mutation'][i]
        ind = 'NA'
        for n_i in name_index:        
            if n_i[0].upper() == pdbid and n_i[1].upper() == mutation:
                ind = n_i[2]
                break
            
        predictable = 0
        for cl in clean_results:
            if cl[3] == ind:
                scores_bASA['comp'][i]=cl[2]
                predictable = 1
                break
        if predictable == 0:
            scores_bASA['comp'][i]='Unpredictable'
    # Save the results
    os.chdir('/home/leo/Documents/Database/Pipeline_New/Codes/Results/Affinity_prediction')    
    with open('clean_results', 'w') as f:
        json.dump(clean_results, f)
    save = scores_bASA[['#pdb', 'mutation', 'expt', 'comp']]  
    save.to_csv('CM.dat', index = False)
#Write_pred_results_in_dat()
#os.chdir('/home/leo/Documents/Database/Pipeline_New/Codes/Results/Affinity_prediction')    
#with open('clean_results', 'r') as f:
#    clean_results = json.load(f)
#len(clean_results)
#####################################################################################
#os.chdir('/home/leo/Documents/Database/Pipeline_New/Results')
#with open('cross_results', 'r') as f:
#    cross_results = json.load(f)
#with open('cross_numerical_Markov_1', 'r') as f:
#    cross_numerical_Markov_1 = json.load(f)
#from pyexcel_ods import get_data, save_data
#data = {'cross_validation_results':cross_results}
#save_data('cross_validation_results.ods', data)
'''
Explanation about the results:
***************************************************
affinity_results are the total results of my prediction
*******************************************************
other_pred_dict is a dictionary of the predictions by other methods on the same set
            of data, it contains the auc, fpr, tpr and CI_95(confidence interval generated by 
            bootstraping)
********************************************************
my_pred_dict is a dictionary of the prediction results by my method. It contains the same
            thing as the other_pred_dict
*******************************************************
auc_all: a list of the auc given by all the methods.
'''
