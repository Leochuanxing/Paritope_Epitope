import copy
import numpy as np
import json
import os
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
        chain_pos = 2
    elif Ab_Ag == 'Ag':
        chain_pos = 3
    for nple in consecutive:
        n_contact = 0
        for i in selected_contact:
            if i[chain_pos] in nple:
                n_contact += i[3]
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
        mutation_match_parameter['mutations']:
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
def Select_contact_opposite(mutation_match_parameter, sequence, cutoff=5):
    # take out the parameter
    pdbid = mutation_match_parameter['pdbid']
    chain = mutation_match_parameter['mutation_chain']
    mutations = mutation_match_parameter['mutations']
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
        if chain == i[0][chain_pos] and i[aa_pos] in mutations \
        and i[0][opposite_chain_pos] == opposite_chain: 
               
            selected_contact.append(i)
            possible_opposite.append(i[opposite_aa_pos]) 
            equal_mutations.append(i[aa_pos])
                
    # We have to make sure the selected_contact contains all the mutations, otherwise our
    # prediction doesn't make sense
    for mut in mutations:
        if mut not in equal_mutations:
            selected_contact = []; possible_opposite = []            
            break        
                             
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
def Paire_select(mutation_match_parameter, sequence):

    # Extract the required information from the pdb file
    selected_contact, possible_opposite, pdbid_sequence = Select_contact_opposite(mutation_match_parameter, sequence)
    # Take out the values from the parameters
    chain = mutation_match_parameter['mutation_chain']
    mutations = mutation_match_parameter['mutations']
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
        for length in [4, 3, 2, 1]:
            longest_possible_consecutives = Get_consecutive(possible_opposite, length, free_type=1)
            if longest_possible_consecutives != []:
                break
            
        for choosen_opposite_pos in longest_possible_consecutives:
            # Correct the direction of the choosen_opposite_pos
            choosen_opposite_pos.sort()
            # define a small function to change the order of the paires
    #        if len(mutations) >= 2 and len(choosen_opposite_pos)>=2:
            if len(choosen_opposite_pos)>=2:
                if Ab_Ag == 'Ab':
                     Order_Ab_sequence(mutations, choosen_opposite_pos, selected_contact)
                elif  Ab_Ag == 'Ag':
                    Order_Ab_sequence(choosen_opposite_pos, mutations, selected_contact)            
    
            # Load the amino acids to the paires according to the choosen_epitope_pos          
            original_aa = []
            for i in mutations:
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
    
    # Load the results
    mutation_match_parameter['paires'] = paires
            
#    return mutation_match_parameter
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
    mutation_aa = mutation_match_parameter['mutations_aa']
    Ab_Ag = mutation_match_parameter['Ab_Ag']
    paires = mutation_match_parameter['paires']
    # Do the work
    if paires != []:
        mutated_paires = []
        if Ab_Ag == 'Ab':
            for parepi in paires:
                mutated_paires.append([mutation_aa, parepi[1]])
        elif Ab_Ag == 'Ag':
            for parepi in paires:
                mutated_paires.append([parepi[0], mutation_aa])
        # Load the results to original_mutation_sets.
        original_mutation_sets = [paires, mutated_paires]
        mutation_match_parameter['original_mutation_sets'] = original_mutation_sets
    else:
        mutation_match_parameter['original_mutation_sets'] = []
    
#    return mutation_match_parameter
           
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
def Predition_RBFN_coverage(testing_set):
#testing_set = testing5
    key = str(len(testing_set[0][0]))+'_'+str(len(testing_set[0][1]))+'_1_1_1_2_1perchain'
    os.chdir('/home/leo/Documents/Database/Pipeline_New/Results')
    with open('test_coverage_RBFN_results', 'r') as f:
        test_coverage_RBFN_results = json.load(f)            

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

def Compare(mutation_match_parameter, model = 'RBFN'):
    # Take out the values
#    fold_change = mutation_match_parameter['fold_change']
    original_mutation_sets = mutation_match_parameter['original_mutation_sets']
    # Get the original sets and the mutated sets
    if original_mutation_sets != []:
        original_sets = original_mutation_sets[0]
        mutated_sets = original_mutation_sets[1]
        # Make predictions for each set
        if model == 'RBFN':
            original_pred = Predition_RBFN_coverage(original_sets)
            mutated_pred = Predition_RBFN_coverage(mutated_sets)
#        if model == 'Logistic':
#            original_pred = Prediction_LogisticRegression(original_sets, wd_results, wd_negative_samples)
#            mutated_pred = Prediction_LogisticRegression(mutated_sets, wd_results, wd_negative_samples)
        original_mutation_score = [np.sum(original_pred), np.sum(mutated_pred)]
    else:
        return 'Empty Match'        
        
    return original_mutation_score

######################################################################
                
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

def Mutation_list_prediction(one_list, sequence,  model = 'RBFN'):  

    # Select the the paires and Generate the original and mutation sets
#    list_parameters = []
    for mutation_match_parameter in one_list:
#        mutation_match_parameter = Paire_select(mutation_match_parameter,sequence)
#        list_parameters.append(copy.deepcopy(Original_mutation_sets(mutation_match_parameter)))
        Paire_select(mutation_match_parameter,sequence)
        Original_mutation_sets(mutation_match_parameter)
    # Make prediction
    # When we make prediction, we simply assume each paratope-epitope paire contribute equally to 
    # the final affinity
    all_original_mutation_score= []
    for parameter in one_list:
        original_mutation_score = Compare( parameter, model)
        # Everytime after Compare we have to switch the working directory back. THis is not good.
        os.chdir('/home/leo/Documents/Database/Pipeline/Complexes with Affinity/structure')
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
    for original_mutaion in all_original_mutation_score:
        scores_original += original_mutaion[0]
        scores_mutation += original_mutaion[1]
    
    return [scores_original, scores_mutation]

######################################################################
#os.chdir('/home/leo/Documents/Database/Pipeline/Complexes with Affinity')
#with open('good_matched_ids', 'r') as f:
#    good_matched_ids = json.load(f)
#with open('good_combined_ids', 'r') as f:
#    good_combined_ids = json.load(f)
#with open('contact', 'r') as f:
#    contact = json.load(f)    
#with open('sequence', 'r') as f:
#    sequence = json.load(f)

#os.chdir('/home/leo/Documents/Database/Pipeline/Complexes with Affinity/structure')



##################################################################
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
def Initiate_mutation_match_parameter(good_matched_ids, good_combined_ids, one_mutation):
    # Take out the values from the one_mutations
    mutations_info = one_mutation['mutations']
    affinities = one_mutation['affinities']
    mutate_to = affinities[0][1]
    # Calculate the fold change
    fold = float(affinities[1][0])/float(affinities[1][1])
    #Load up the fold_change
    fold_change = [affinities[1][0], affinities[1][1], affinities[2], fold]
    # Finde the pdbid
    pdbid = mutations_info[0][0]
    # Find the combined ids
    combined_ids = good_combined_ids[pdbid]
    # Find the matched ids
    matched_ids = good_matched_ids[pdbid]
    # Get the opposite chain
    # create an empty list
    unit_list = []
    # Load up the list
    for i in range(len(mutations_info)):
        sub_mutaion = mutations_info[i]
        # Find the mutaion chain
        mutation_chain = sub_mutaion[1]    
        # Find the value of Ab_Ag and opposite_chains
        opposite_chains = []
        for match in matched_ids:            
            for i in range(3):
                if match[i] == mutation_chain and i == 2:
                    Ab_Ag = 'Ag'
                    opposite_chains.extend([match[0], match[1]])                  
                elif match[i] == mutation_chain and i != 2:
                    Ab_Ag = 'Ab'
                    opposite_chains.append(match[2])
        # Find the opposite chains
#        opposite_chains = sub_mutaion[4]
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
            temp_mutation_match_parameter['mutations'] = sub_mutaion[2]
            # Load the mutations_aa
            temp_mutation_match_parameter['mutations_aa'] = sub_mutaion[3]
            # Load the opposite chain 
            temp_mutation_match_parameter['opposite_chain'] = opposite
            # Load the mutate to
            temp_mutation_match_parameter['mutate_to'] = mutate_to
            # Load the fold change
            temp_mutation_match_parameter['fold_change'] = fold_change
            # Load the matched_ids
            temp_mutation_match_parameter['matched_ids'] = matched_ids
            # Load the combined_ids
            temp_mutation_match_parameter['combined_ids'] = combined_ids            
            
            # Add the dictionary to the unit_list
            unit_list.append(copy.deepcopy(temp_mutation_match_parameter))
    
    return unit_list
############################################################
'''
List_list_other:
    a function to generate similar output as List_list, the only difference is that 
    the input is a little bit different from List_list
Input:
    other_mutations: mutations collected from papers
    good_matched_ids:
    good_combined_ids:
    
Output:
    list_list_other:
        in the same form as the output of List_list
'''
def List_list_other(other_mutations, good_matched_ids, good_combined_ids):
    list_list_other = []
    # Initiate the mutation_match_parameters
    for one_mutation in other_mutations:
        unit_list = Initiate_mutation_match_parameter(good_matched_ids, good_combined_ids, one_mutation)
#        for mutation_match_parameter in unit_list:
#            Paire_select(mutation_match_parameter)
#            Original_mutation_sets(mutation_match_parameter)
#            # Load to the list_list_other
        list_list_other.append(copy.deepcopy(unit_list))
        
    return list_list_other
    
#################################################################

'''
Predict_list_list:
    a function to make prediction on the given list_list
Input:
    list_list:
        a list of  basic prediction unit, which is list of mutation_match_parameter dictionaries
Output:
    results_list: a list, each element is a 
        a dictionary in the following form:
        results['original'] = original pdbid
        results['mutate_to'] = the pidbid of the complex the original complex is mutated to
        results['fold_change'] = the fold change of the affinity  original_affinity/mutation_affinity
        results['predictions'] = [original_prediction_value, mutation_prediction_value]
'''
def Predict_list_list(list_list, wd_results, wd_negative_samples,sequence, model):
    results_list = []
    for i in range(len(list_list)):
        list_1 = list_list[i]
        # Change the working directory
        os.chdir('/home/leo/Documents/Database/Pipeline/All with peptide 5+ resolution 4A/structure')
               
        # Create a dictionary to contain the detailed prediction results
        results_1 = {}
        # Make prediction for the lists
        if list_1 != []:
            prediction_1 = Mutation_list_prediction(list_1, wd_results, wd_negative_samples,sequence, model)
            # Store the results
            results_1['original'] = list_1[0]['pdbid']
            results_1['mutate_to'] = list_1[0]['mutate_to']
            results_1['fold_change'] = list_1[0]['fold_change']
            results_1['predictions'] = prediction_1
       
        if results_1 != {}:
            results_list.append(results_1)
    return results_list

#####################################################################
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
    if affinity_type == 'Kd':
        if fold > 1 and prediction1 < prediction2:
            right_or_wrong = True
        elif fold < 1 and prediction1 > prediction2:
            right_or_wrong = True
        else:
            right_or_wrong = False
            
    # In the case when the affinity_type is Ka      
    elif affinity_type == 'Ka':
        if fold > 1 and prediction1 > prediction2:
            right_or_wrong = True
        elif fold < 1 and prediction1 < prediction2:
            right_or_wrong = True
        else:
            right_or_wrong = False
    
    return right_or_wrong
########################################################################
'''
Count the number of predictions and the number of correct predictions
'''
def Count_correct(results_list, fold_change_cut = 1):
    total_prediction = 0
    correct_prediction = 0
    
    for results in results_list:
        # Take out the values
        fold_change = results['fold_change']
        fold = fold_change[3]
        affinity_type = fold_change[2]
        predictions = results['predictions']
        
        if predictions[0] != 'Empty Match' and fold >= fold_change_cut: 
            total_prediction += 1
            if Right_or_wrong(fold, affinity_type, predictions[0], predictions[1]):
                correct_prediction += 1
                
        elif predictions[0] != 'Empty Match' and fold <= 1/fold_change_cut:
            total_prediction += 1
            if Right_or_wrong(fold, affinity_type, predictions[0], predictions[1]):
                correct_prediction += 1

    return total_prediction, correct_prediction
###############################################################################
'''
Do the prediction, let us use one paires of complexes as an example to build up the 
steps
'''
#  Assign the working directories
directory_paires = [['/home/leo/Documents/Database/Pipeline/Results/1_free',\
  '/home/leo/Documents/Database/Pipeline/Negative samples/1_free'],\
 ['/home/leo/Documents/Database/Pipeline/Results/0_free',\
                       '/home/leo/Documents/Database/Pipeline/Negative samples/0_free']]
        
wd_results = directory_paires[0][0]
wd_negative_samples = directory_paires[0][1]

# Open the workable list
os.chdir('/home/leo/Documents/Database/Pipeline/Affinity/All_structures')
with open('mutations', 'r') as f:
    mutations = json.load(f)
with open('combined_ids', 'r') as f:
    good_combined_ids = json.load(f)
with open('matched_ids', 'r') as f:
    good_matched_ids = json.load(f)
with open('mutations_mix', 'r') as f:
    mutations_mix = json.load(f)
with open('sequence', 'r') as f:
    sequence = json.load(f)
with open('mutations_1dvf', 'r') as f:
    mutations_1dvf = json.load(f)
#len(mutations)
############################################################
if __name__ == '__main__':
    list_list = List_list_other(mutations, good_matched_ids, good_combined_ids)
    results_list = Predict_list_list(list_list, wd_results, wd_negative_samples,sequence, model = 'RBFN')
    total_prediction, correct_prediction = Count_correct(results_list, fold_change_cut = 1)   
    total_prediction
    correct_prediction
    results_list
    
total_prediction
correct_prediction
40/49

