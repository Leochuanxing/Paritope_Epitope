'''PREDICT MUTATION'''
import json
#import xlrd
import os
import copy
import numpy as np
#from pyexcel_ods import get_data, save_data
from Extract_mutation_pairs import Formalized_contacting
'''#############################################################'''
from Bio import Align
from Bio.SubsMat.MatrixInfo import blosum62
aligner = Align.PairwiseAligner()
aligner.substitution_matrix = blosum62
aligner.mode = 'global'
'''###################################################################'''
os.chdir('/home/leo/Documents/Database/Data_Code_Publish/Codes')
from AlignmentConstraint import To_seq
from RBFN_coverage import Distance_matrix, Design_matrix, Observed_values
from RBFN_coverage import  Coverage_reduce_centers, NumpyEncoder
'''######################################################################'''
# First change the mutation into the inputs of function Complete_extract_contacting
# The basic strategy is to process each mutaion set with measured affinity as a unit
# and iterate all the mutaion set
'''  
search_para = {}          
search_para['moving'] 
search_para['step_size']
search_para['start_dist']
search_para['end_dist']
search_para['cut_dist']
search_para['form']
search_para['within_range']
search_para['pdbid']
earch_para['mut_pos']
search_para['mut_chain_id']
'''

'''##############################################################
Processing_mutation_set: a sub function of Workable input
'''

def Processing_mutation_set(mutation_set):
    # Separate the information
    DDG = mutation_set[-1][1]
    pdbid = mutation_set[0][0]
    mut_id = mutation_set[0][5]
    mutations = mutation_set[:-1]
    
    # Get the chains, and pos
    chains = []
    for mu in mutations:
        if mu[1] not in chains:
            chains.append(mu[1])
    
    mut_info = []
    for chain in chains:
        pos = [mu[2] for mu in mutations if mu[1] == chain]
        mut_aa = [mu[3] for mu in mutations if mu[1] == chain]

        mut_info.append(copy.deepcopy([pdbid, chain, pos, mut_aa, DDG, mut_id]))
    
    return mut_info
# Read the parsed mutation information
'''
Output:
    workable_input:
        a list with elements in the form 
        
        [...,
        [pdbid, mut_chain_id, [pos], [mut_aa], affinities, mut_id],
        ...]
        
        The element is a list with the number of elements equals the number of 
        mutation chains.
'''
def Workable_input(mutation_d):
    os.chdir(mutation_d)
    formated = get_data('Formated.ods')

    workable_input = []
    for pdb, mutations in formated.items():
        begin = 0; end = 0
        for ind, mut in enumerate(mutations):
            if mut[0] == 'affinities':
                end = ind + 1
                mutation_set = mutations[begin: end]
                begin = end
                # processing the mutation set
                mut_info = Processing_mutation_set(mutation_set)
                workable_input.append(copy.deepcopy(mut_info))
            
    return workable_input
  
'''##############################################################3'''
'''
Generate_wt_mut_one: a function to generat wt_mut pairs for one_tuple
Input:
    one_tuple:
        in the form
        
        ( mut_pos, mut_aa, mut_chain_id,  op_pos, op_aa, op_chain_id,\
        contact_number, mut_chain_type, op_chain_type)
        
        if form == 'flanked', mut_aa and op_aa are lists
Output:
    Ab_Ag_wt_mut_pair: a list in the form
    [[wt_Ab_Ag_pair], [mut_Ab_Ag_pair]]
'''
def Generate_wt_mut_one(one_tuple, form, mut_aa):
    mut_chain_type = one_tuple[7]
    # Extract wt information from one_tuple and generate mut_info
    if form == 'one' or 'multiple':
        wt_pair = [[one_tuple[1]], [one_tuple[4]]]
        mut_pair = [[mut_aa], [one_tuple[4]]]

    elif form == 'flanked':
        wt_pair = [one_tuple[1], one_tuple[4]]
        
        tri = one_tuple[1]
        tri[1] = mut_aa
        mut_pair = [tri, one_tuple[4]]
        
    # Put the Ab first
    if mut_chain_type == 'A':
        Ab_Ag_wt_pair = [wt_pair[1], wt_pair[0]]
        Ab_Ag_mut_pair = [mut_pair[1], mut_pair[0]]
    else:
        Ab_Ag_wt_pair = wt_pair
        Ab_Ag_mut_pair = mut_pair
        
    
    Ab_Ag_wt_mut_pair = [Ab_Ag_wt_pair, Ab_Ag_mut_pair]
    
    return Ab_Ag_wt_mut_pair
'''
One_chian_output: a sub fuction of Workable_output
Output:
    a list with elements in the following form: 
    
    [[wt_Ab_Ag_pair], [mut_Ab_Ag_pair],contact_number, DDG, mut_id]
    
'''
def One_chian_output(one_chain, search_para, combined_ids, sequence, structure_d):
    
    search_para['pdbid'] = one_chain[0]
    search_para['mut_chain_id'] = one_chain[1]
    search_para['mut_pos'] = one_chain[2]
    
    mut_pos_list = one_chain[2]
    mut_aa_list = one_chain[3]
    DDG = one_chain[4]
    mut_id = one_chain[5]
    
    form = search_para['form']
    
    formalized_pairs = \
        Formalized_contacting(search_para, combined_ids, sequence, structure_d)
    
    # reformat the formalized pairs with more information and transform into wt mut pairs
    Ab_Ag_wt_mut_pairs_list = []
    for one_tuple in formalized_pairs:
        for pos, aa in zip(mut_pos_list, mut_aa_list):
            if one_tuple[0] == pos:
                mut_aa = aa; contact_number = one_tuple[6]
                Ab_Ag_wt_mut_pair = Generate_wt_mut_one(one_tuple, form, mut_aa)
                # Attach the DDG and the mutid
                Ab_Ag_wt_mut_pair.extend([contact_number, DDG, mut_id])
                # Load to the Ab_Ag_wt_mut_pairs_list
                Ab_Ag_wt_mut_pairs_list.append(copy.deepcopy(Ab_Ag_wt_mut_pair))

        
    return Ab_Ag_wt_mut_pairs_list
'''
Input:
    search parameter contains the following values:
        search_para['moving'] 
        search_para['step_size']
        search_para['start_dist']
        search_para['end_dist']
        search_para['cut_dist']
        search_para['form']
        search_para['within_range']
        
Output:
    workable, a list with elements lists in the following form
    
    [...'
    [[wt_Ab_Ag_pair], [mut_Ab_Ag_pair], contact_number, DDG,mut_id],
    ...]
    

'''
def Workable_output(workable_input, search_para, combined_ids, sequence, structure_d):
    workable = []

    for one_mutation in workable_input:
        one_mut_workable = []
        for one_chain in one_mutation:
            Ab_Ag_wt_mut_pairs_list =\
                One_chian_output(one_chain, search_para, combined_ids, sequence, structure_d)
            one_mut_workable.extend(Ab_Ag_wt_mut_pairs_list)
            
        workable.append(copy.deepcopy(one_mut_workable))
    
    return workable

'''############################################################################'''
# Give a little test of the workable output


'''#################################################################################'''
def Sub_predict(wt_pair, mut_pair, model_all):
    
    # Calculate the design matrix    
    prefix = str(len(wt_pair[0])) + '_' + str(len(wt_pair[1]))
    key = prefix + '_0_0_1_2_1perchain'
    # Get the centers and coefficients for related match type
    centers = model_all[key]['centers']
    coefficients = model_all[key]['coefficients']
    
    # Calculate the design matrix
    distance_matrix = Distance_matrix([wt_pair, mut_pair], centers, square=False)
    radius_coeff = np.ones((len(centers), 1))
    design_matrix = Design_matrix(distance_matrix, radius_coeff)
    # Attach the constant column
    design_matrix = np.hstack((design_matrix, np.ones((2, 1))))
    
    # Predict
    coefficients = np.reshape(coefficients, (-1, 1))
    pred = design_matrix.dot(coefficients)
    
    # Take the difference between the two predicted values
    diff = pred[1] - pred[0]
    
    return diff

###############################################################################
'''
working_d: Where the train_results locates.
'''



def Predict_affinity(workable, working_d, binary = True):
    predict_affinity_results = []
    # Open the training results
    os.chdir(working_d)
    with open('train_results', 'r') as f:
        train_results = json.load(f)
        
    # Set the model
    if binary:
        model_all = train_results['cross_binary_Gaussian_']
        best_gap_ext = model_all['best_gap_ext']
        aligner.gap_score = best_gap_ext[0]
        aligner.extend_gap_score = best_gap_ext[1]
    else:
        model_all = train_results['cross_numerical_Gaussian_']
        best_gap_ext = model_all['best_gap_ext']
        aligner.gap_score = best_gap_ext[0]
        aligner.extend_gap_score = best_gap_ext[1]
        
    for one_mut_set in workable:

        if one_mut_set != []:
            # Calculate the sum of the contact number
            contact_sum = 0; weighted_diff = 0
            for one_pair in one_mut_set:
                wt_pair = one_pair[0]
                mut_pair = one_pair[1]
                diff = Sub_predict(wt_pair, mut_pair, model_all)
                # Take the weighted average
                contact_sum += one_pair[2]
                DDG = one_pair[3]
                mut_id = one_pair[4]

                
                weighted_diff += diff * one_pair[3]
                
            weighted_diff /= contact_sum
            
            # Load to the predict_affinity_results
            predict_affinity_results.append([mut_id, DDG, weighted_diff[0]])
        else:
            predict_affinity_results.append(['', '', 'Unpredicable'])
    
    return predict_affinity_results


'''##########################################################################'''
def Calculate_AUC(pred, observed_values):
    # Calculate the number of positive observations
    positive_total = 0 ; negative_total = 0
    for observed in observed_values:
        if observed > 0 :
            positive_total += 1
        elif observed < 0:
            negative_total += 1
    if negative_total == 0:
        print('None negative samples in Calculating AUC')
        return None, None, None
    if positive_total == 0:
        print('None positive samples in Calculating AUC')
        return None, None, None
    # Match_up the pred with the observed_values   
    match_up = []
    for j in range(len(observed_values)):
        match_up.append([pred[j], observed_values[j]])
    match_up.sort(key = lambda x:x[0], reverse = True)

    # Calculate the FPR and TPR 
    FPR = []; TPR = []
    n_positive = 0
    n_negative = 0
    for match in match_up:
        if match[1] > 0:
            '''Here we use >0 to make sure it works for both binary and non binary'''
            n_positive += 1
        elif match[1] < 0:
            n_negative += 1
            
        TPR.append(n_positive/positive_total)
        FPR.append(n_negative/negative_total)
    
    # Calculate the AUC
    AUC = 0
    for i in range(1, len(TPR)):
        AUC += TPR[i] * (FPR[i] - FPR[i-1])
    
    return AUC, TPR, FPR
#################################################################################
'''
Input:
    cut_DDG_lower, cut_DDG_upper: cut values, only to consider the mutations with 
            the DDG in 
                        (cut_DDG_lower, cut_DDG_upper]
            open in the left and closed in the right.
'''
def Analyze_resutls(predict_affinity_results, cut_DDG_lower, cut_DDG_upper):
    # Select the prediction by cut_DDG
    selected_cut_DDG = []
    for res in predict_affinity_results:
        if res[2] != 'Unpredicable' and res[1] != 0:
            if abs(res[1]) > cut_DDG_lower and abs(res[1]) <= cut_DDG_upper:
                selected_cut_DDG.append(res)
    # Calculate AUC, FPR, TPR
    observed_values = []; pred = []
    for match in selected_cut_DDG:
        # We have to add a nagetive sign to the observed to make it consistant 
        observed_values.append(-match[1])
        pred.append(match[2])
        
    AUC, TPR, FPR = Calculate_AUC(pred, observed_values)
    
    # Calculate the absolute correct
    number_correct = 0
    for match in selected_cut_DDG:
        if match[1] * match[2] < 0:
            number_correct += 1
    correct_ratio = number_correct / len(selected_cut_DDG)
    
    return selected_cut_DDG, AUC, TPR, FPR, correct_ratio

'''####################################################################'''

if __name__ == '__main__':  

    preliminary_pred = {}
      
    mutation_d = '/home/leo/Documents/Database/Data_Code_Publish/Mutations'
    workable_input = Workable_input(mutation_d) 
    
    os.chdir('/home/leo/Documents/Database/Data_Code_Publish/Structures')
    with open('sequence', 'r') as f:
        sequence = json.load(f)
    with open('combined_ids', 'r') as f:
        combined_ids = json.load(f)
        
    search_para = {}
    search_para['moving'] = True
    search_para['step_size'] = 0.25
    search_para['start_dist'] = 3
    search_para['end_dist'] = 8
    search_para['cut_dist'] = 6
    
    for form in ['one', 'multiple']:        
        search_para['form'] = form
        for within_range in [True, False]:
            search_para['within_range'] = within_range
            
            print('Working on: '+ form +'    '+ str(within_range))
            
            # Container
            container = []

            structure_d = '/home/leo/Documents/Database/Data_Code_Publish/Structures/imgt'    
            workable = Workable_output(workable_input, search_para, combined_ids, sequence, structure_d)
            
            working_d = '/home/leo/Documents/Database/Data_Code_Publish/Codes/Results'
            predict_affinity_results = Predict_affinity(workable, working_d, binary= True)
            
            for ran in [[0,0.5], [0, 100], [0.5, 100], [1, 100]]:
                cut_DDG_lower = ran[0]
                cut_DDG_upper = ran[1]
                selected_cut_DDG, AUC, TPR, FPR, correct_ratio = \
                            Analyze_resutls(predict_affinity_results, cut_DDG_lower, cut_DDG_upper)
                container.append([ran, AUC, len(selected_cut_DDG)])
                
            preliminary_pred[form+'_WithinRange_'+str(within_range)] = copy.deepcopy(container)
            
            saving_d = '/home/leo/Documents/Database/Data_Code_Publish/Codes/Results'
            os.chdir(saving_d)
            with open('affinity_pre_results', 'w') as f:
                json.dump(preliminary_pred, f)
                

#correct_ratio
#AUC
#len(TPR)
#len(selected_cut_DDG)
#selected_cut_DDG

#os.chdir('/home/leo/Documents/Database/Data_Code_Publish/Codes/Results')
#with open('train_results', 'r') as f:
#    train_results = json.load(f)

#train_results.keys()
#train_results['cross_binary_Gaussian_'].keys()
#train_results['cross_binary_Gaussian_']['1_1_0_0_1_2_1perchain'].keys()
#len(train_results['cross_binary_Gaussian_']['1_1_0_0_1_2_1perchain']['centers'])
#len(train_results['cross_binary_Gaussian_']['1_1_0_0_1_2_1perchain']['coefficients'])
#
#train_results['cross_binary_Gaussian_']['1_1_0_0_1_2_1perchain']['centers'][:5]


























