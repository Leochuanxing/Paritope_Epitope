import copy
import numpy as np
import json
import os
from pyexcel_ods import get_data, save_data
#import math
os.chdir('/home/leo/Documents/Database/Pipeline_New/Codes')
from AAC_2 import Coordinates, Get_contact
from FrameConstraint import Get_consecutive
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

#########################################################################################

'''
Criteria:
    First we pick the number of mutations positions, they should be continuous in the 1_free sense.
    
    Second, we find all the opposite consecutive antigen amino acids in the 1_free sense.
             Among all those consecutive antigen amino acids with length no more than 4, we 
             choose the one with the largest contact numbers with the given mutation positions.   
'''
###################################################################################
      
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
def Select_contact_opposite(mutation_match_parameter, sequence, cutoff=7):
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
        cdn = Coordinates(file, combined_ids)
#    with open(pdbid+'.pdb', 'r') as file:
#        pdbid_sequence = Chain_seq(file, combined_ids)
    pdbid_sequence = sequence[pdbid]
    # store the paires in paires
    possible_opposite = []
    if Ab_Ag == 'Ab':
        chain_pos = 2; opposite_chain_pos = 3; aa_pos = 1; opposite_aa_pos = 2
    elif Ab_Ag == 'Ag':
        chain_pos = 3; opposite_chain_pos = 2; aa_pos = 2; opposite_aa_pos = 1
        
    contact = Get_contact(cdn, matched_ids, cutoff)
    # Carry out the above process:
    # take out all the contact containing the chain_name
    selected_contact = []; possible_opposite = []; equal_mutations = []    
    for i in contact:
        if chain == i[0][chain_pos] and i[aa_pos] in mutations_pos \
        and i[0][opposite_chain_pos] == opposite_chain: 
               
            selected_contact.append(i)
            possible_opposite.append(i[opposite_aa_pos]) 
            equal_mutations.append(i[aa_pos])
                
    # We have to make sure the selected_contact contains all the mutations, otherwise our
    # prediction doesn't make sense
    '''*********************************'''
    for mut in mutations_pos:
        if mut not in equal_mutations:
            selected_contact = []; possible_opposite = []            
            break        
       
    '''*******************************'''                      
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
def Paire_select(mutation_match_parameter, sequence, mode = 'single', cutoff = 6):

    # Extract the required information from the pdb file
    selected_contact, possible_opposite, pdbid_sequence =\
        Select_contact_opposite(mutation_match_parameter, sequence, cutoff)
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
        for length in [3,2,1]:
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
def Predition_RBFN_coverage(test_coverage_RBFN_results, testing_set):
#testing_set = testing5
    key = str(len(testing_set[0][0]))+'_'+str(len(testing_set[0][1]))+'_0_0_1_2_1perchain'    

    non_redundent_training_set = test_coverage_RBFN_results[key]['non_redundent_training_set']
    distance_matrix = Distance_matrix(testing_set, non_redundent_training_set, square = False)
    testing_design_matrix = Design_matrix(distance_matrix)
    
    coeff = test_coverage_RBFN_results[key]['coeff']
    coeff = np.array(coeff).reshape((-1,1))
    #Calculate the prediction results
    predictions = testing_design_matrix.dot(coeff)
    
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
def Sub_parse_mutations(mutations_list_onepdb, free_type):
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
def Parse_mutations(mutation_data, free_type):
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
        # get the list of mutaions
        mutations_list_onepdb = mutation_data[start:end]    
        # grouped mutations by free_type
        mutations_by_free_type_chain = Sub_parse_mutations(mutations_list_onepdb, free_type)
        # affinities
        affinities = mutation_data[end]
        
        separated_mutation_complex.append([pdbid, mutations_by_free_type_chain,affinities])
        
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
        a dictionary, with keys 'mutations', 'affinities'. This dictionay is one of
        the element of the recorded mutaions from the literatures
Output:
    unit_list:
        a list of mutation_match_parameters, and this list is a basic prediction unit
'''

def Initiate_mutation_match_parameter(matched_ids, combined_ids, one_mutation):
    # Take out the values from the one_mutations
    pdbid = one_mutation[0]    
    # Find the affinity type    
    affinity_type = one_mutation[-1][3]    
    # Calculate the fold change
    '''Here we would like to transform all the fold change to Kd '''
    if affinity_type == 'Kd':
        fold_change = one_mutation[-1][2]/one_mutation[-1][1]
    elif affinity_type == 'Ka':
        fold_change = one_mutation[-1][1]/one_mutation[-1][2]
    elif affinity_type == 'DDG':
        if one_mutation[-1][2] == 0:
            fold_change = 1
        else:
            fold_change = np.exp(one_mutation[-1][2]/(8.31*298))
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

def Mutation_list_prediction(one_list, sequence,  test_coverage_RBFN_results, mode = 'single', cutoff=6):  

    # Select the the paires and Generate the original and mutation sets
#    list_parameters = []
    for mutation_match_parameter in one_list:
        Paire_select(mutation_match_parameter,sequence, mode, cutoff)
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
        scores_original += original_mutaion[0]/length
        scores_mutation += original_mutaion[1]/length
    
    return [scores_original, scores_mutation]

######################################################################
'''
Test the above functions
'''
#data = small_data
def Predict(data, matched_ids, combined_ids, sequence, test_coverage_RBFN_results, mode = 'single', cutoff=7):
    mu_list = Parse_mutations(data, free_type = 0)
    results = []

    for mutation in mu_list: 
        sub_result = []
        one_list = Initiate_mutation_match_parameter(matched_ids, combined_ids, mutation)
        scores = Mutation_list_prediction(one_list, sequence,  test_coverage_RBFN_results, mode, cutoff)
        sub_result.append(mutation[0])
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
            
#    # In the case when the affinity_type is Ka      
#    elif affinity_type == 'Ka':
#        if fold > 1 and prediction1 > prediction2:
#            right_or_wrong = True
#        elif fold < 1 and prediction1 < prediction2:
#            right_or_wrong = True
#        else:
#            right_or_wrong = False
#    elif affinity_type == 'DDG':
#        fold = np.exp(fold*1000/(8.31*298))
#        if fold > 1 and prediction1 > prediction2:
#            right_or_wrong = True
#        elif fold < 1 and prediction1 < prediction2:
#            right_or_wrong = True
#        else:
#            right_or_wrong = False    
    
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
def ROC_AUC(results_all, fold_cut = 1):
    res_all = []
    for results in results_all:
        
        if results[-1][0]!= 'Empty Match' and  results[-1][2] != 1:
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
    tp = 0
    fp = 0
    for res in res_all:
        if res[0] < 1/fold_cut:
            tp += 1
        elif res[0] > fold_cut:
            fp += 1
        TPR.append(tp/total_positive)
        FPR.append(fp/total_negative)
    # Calculate AUC
    AUC = 0
    for i in range(1, len(TPR)):
        AUC += TPR[i] * (FPR[i]- FPR[i-1])
    
    return TPR, FPR, AUC
    
    
################################################################################
'''
[['1vfb', 'B', 27, 'ASN'],
 ['affinities', 1.49, 1.9, 'Kd'],
 [],
 ['1vfb', 'B', 27, 'GLN'],
 ['affinities', 1.49, 1.7, 'Kd']]
'''
# Change the mutations into the forms given above   

###########################################################################

if __name__ == '__main__':
    os.chdir('/home/leo/Documents/Database/Pipeline_New/All with peptide 5+ resolution 4A')
    with open('matched_ids', 'r') as f:
        matched_ids = json.load(f)
    with open('combined_ids', 'r') as f:
        combined_ids = json.load(f)
    with open('sequence', 'r') as f:
        sequence = json.load(f)
        
    os.chdir('/home/leo/Documents/Database/Pipeline_New/Results')
    with open('test_coverage_RBFN_0_0_results', 'r') as f:
        test_coverage_RBFN_results = json.load(f)
        
    os.chdir('/home/leo/Documents/Database/Pipeline_New/Results')
    data_raw = get_data('Affinities.ods')
    keys = data_raw.keys()
#    keys = ['1jrh']
    affinity_results_dict_single = {}
    for cut in np.arange(4, 8.5, 0.5):
        results = []
        for key in keys:
            print('working on:  ' + key)
            data = []
            for d in data_raw[key]:
                if d != []:
                    data.append(d)
    
            os.chdir('/home/leo/Documents/Database/Pipeline_New/All with peptide 5+ resolution 4A/structure')
            
            results.extend(Predict(data, matched_ids, combined_ids, sequence,\
                              test_coverage_RBFN_results, mode = 'single', cutoff=cut))
            
        affinity_results_dict_single[str(cut)] = copy.deepcopy(results)
###############################################################################
class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)
###########################################################
os.chdir('/home/leo/Documents/Database/Pipeline_New/Codes/Results')
with open('affinity_results_dict_single', 'w') as f:
    json.dump(affinity_results_dict_single, f, cls = NumpyEncoder)

#    TPR, FPR, AUC = ROC_AUC(results)

#TPR, FPR, AUC = ROC_AUC(results, fold_cut=2.5)

#len(results)
#TPR, FPR, AUC = ROC_AUC(results)     
#AUC  
#        total_prediction, correct_prediction = \
#                        Count_correct(results, fold_change_cut = 1)
    
    # Save the results
#    os.chdir('/home/leo/Documents/Database/Pipeline_New/Results')
#    with open('affinity_testing_results', 'w') as f:
#        json.dump(results, f)

#print(len(data), total_prediction, correct_prediction)

#results

#TPR, FPR, AUC = ROC_AUC(results)

#os.chdir('/home/leo/Documents/Database/Pipeline_New/Results')
#data_raw = get_data('Affinities.ods')
#keys = data_raw.keys()
#for key in keys:
#    data = data_raw[key]
#    for mu in data:
#        if mu != [] and mu[3] != '':
#            if mu[3] == 'DDG':
#                mu[2] = mu[1]
#                mu[1] = 0
#        data_raw.update({key: data})
#save_data('Affinities.ods', data_raw)
        
#    
#data = []
#for d in data_raw:
#    if d != []:
#        data.append(d)
