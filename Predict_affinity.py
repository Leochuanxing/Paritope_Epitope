'''PREDICT MUTATION'''
import sys
import json
import math
import os
import copy
import random
import numpy as np
import pandas as pd
from pyexcel_ods import get_data, save_data


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
from Extract_mutation_pairs import Formalized_contacting
from Some_basic_fuctions import AUC_TPR_FPR, Bootstrap_AUC
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
    mut_info = []
    # Separate the information
    DDG = mutation_set[-1][1]
    pdbid = mutation_set[0][0]
    mut_id = mutation_set[0][5]
    mutations = mutation_set[:-1]
    
    # Get the chains, and pos
#    chains = []
    for mu in mutations:
        mut_info.append([pdbid, mu[1], [mu[2]], [mu[3]], DDG, mut_id])
#        if mu[1] not in chains:
#            chains.append(mu[1])
#    
#    mut_info = []
#    for chain in chains:
#        pos = [mu[2] for mu in mutations if mu[1] == chain]
#        mut_aa = [mu[3] for mu in mutations if mu[1] == chain]
#
#        mut_info.append([pdbid, chain, pos, mut_aa, DDG, mut_id])
    
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
                workable_input.append(mut_info)
            
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
    if form == 'one' or form == 'multiple':
        wt_pair = [[one_tuple[1]], [one_tuple[4]]]
        mut_pair = [[mut_aa], [one_tuple[4]]]

    elif form == 'flanked':
        wt_pair_1 = [x for x in one_tuple[1] if x != '']
        wt_pair_2 = [x for x in one_tuple[4] if x != '']
        wt_pair = [wt_pair_1 , wt_pair_2]
        
        tri = one_tuple[1][:]
        tri[1] = mut_aa
        mut_pair_1 = [x for x in tri if x != '']
        mut_pair_2 = [x for x in one_tuple[4] if x!= '']
        mut_pair = [mut_pair_1 , mut_pair_2]
        
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
def Single_mut_output(single_mut, search_para, combined_ids, matched_ids, sequence, structure_d):
    
    search_para['pdbid'] = single_mut[0]
    search_para['mut_chain_id'] = single_mut[1]
    search_para['mut_pos'] = single_mut[2]
    
    mut_pos_list = single_mut[2]
    mut_aa_list = single_mut[3]
    DDG = single_mut[4]
    mut_id = single_mut[5]
    
    form = search_para['form']
    
    formalized_pairs = \
        Formalized_contacting(search_para, combined_ids, matched_ids, sequence, structure_d)
    
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
    [[wt_Ab_Ag_pair], [mut_Ab_Ag_pair], contact_number, DDG, mut_id],
    ...]
    

'''
def Workable_output(workable_input, search_para, combined_ids,matched_ids, sequence, structure_d):
    workable = []

    for one_mutation in workable_input:
        one_mut_workable = []
        for single_mut in one_mutation:
            Ab_Ag_wt_mut_pairs_list =\
                Single_mut_output(single_mut, search_para, combined_ids, matched_ids, sequence, structure_d)
            one_mut_workable.extend(Ab_Ag_wt_mut_pairs_list)
            
        workable.append(copy.deepcopy(one_mut_workable))
    
    return workable

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
####################################################################################
'''THIS PREDICTION IS TO USE THE PREVIOUS RESULTS'''
def Predition_RBFN_coverage(wt_pair, mut_pair, test_results):
#testing_set = testing5
    key = str(len(wt_pair[0]))+'_'+str(len(wt_pair[1]))+'_0_0_1_2_1perchain' 
    '''********************Pay attention to this scale factor****************'''
#    scale_factor = len(testing_set[0][0]) * len(testing_set[0][1])

    coeff = test_results[key]['coeff']
    coeff = np.reshape(np.array(coeff), (-1,1))
    centers_selected = test_results[key]['centers_selected']
    # Calculate the distance matrix
    distance_matrix = Distance_matrix([wt_pair, mut_pair], centers_selected, square = False)
    # Calculate the design matrix
    linear_coeff = coeff[:len(centers_selected)+1, 0]
    linear_coeff = np.reshape(linear_coeff, (-1,1))
    radius_coeff = coeff[len(centers_selected)+1:, 0]
    radius_coeff = np.reshape(radius_coeff, (-1, 1))
    
    testing_design_matrix = Design_matrix(distance_matrix, radius_coeff, basis_function = 'Gaussian')
    testing_design_matrix = np.hstack((testing_design_matrix, np.ones((2, 1))))
    
    pred = testing_design_matrix.dot(linear_coeff)
    
    diff = pred[1] - pred[0]
    '''******************************'''
#    predictions *= scale_factor
    
    return diff

###############################################################################
'''
working_d: Where the train_results locates.
Output:
    predict_affinity_results, a list with elements in the following form
    
        [mut_id, DDG, weighted_diff]
    
'''



def Predict_affinity(workable, working_d, binary = True):
    predict_affinity_results = []
    # Open the training results
    os.chdir(working_d)
    with open('train_results', 'r') as f:
        train_results = json.load(f)
#   #######################################################################     
#    os.chdir('/home/leo/Documents/Database/Pipeline_New/Complexes/Results')
#    with open('test_results', 'r') as f:
#        test_results = json.load(f)
#    aligner.gap_score = -5
#    aligner.extend_gap_score = -1
#############################################################################33    
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
##################################################################################        
    for one_mut_set in workable:

        if one_mut_set != []:
            print('Working on   ', one_mut_set[0][-1])
            # Calculate the sum of the contact number
            contact_sum = 0; weighted_diff = 0
            for one_pair in one_mut_set:
                wt_pair = one_pair[0]
                mut_pair = one_pair[1]
                try:
#                    diff = Predition_RBFN_coverage(wt_pair, mut_pair, test_results)
                    diff = Sub_predict(wt_pair, mut_pair, model_all)
                except:
                    print(wt_pair, mut_pair)
                    sys.exit()
                # Take the weighted average
                contact_sum += one_pair[2]
                DDG = one_pair[3]
                mut_id = one_pair[4]
                
#                weighted_diff += diff * one_pair[2]
                weighted_diff += diff
                
#            weighted_diff /= contact_sum
            weighted_diff /= len(one_mut_set)
            
            # Load to the predict_affinity_results
            predict_affinity_results.append([mut_id, DDG, weighted_diff[0]])
        else:
            predict_affinity_results.append(['', '', 'Unpredicable'])
    
    return predict_affinity_results


'''##########################################################################'''
#################################################################################
def Calculate_concentration(TPR, top_percent = 0.1):
    concentration = round(TPR[math.floor(len(TPR) * top_percent)] / top_percent, 2)
    
    return concentration
################################################################################
def Calculate_correlation(selected_cut_DDG):
    ar = np.array(selected_cut_DDG)[:,[1,2]]
    ar_mean = np.average(ar, axis=0).reshape((1,-1))
    ar_d_mean = ar - ar_mean
    x_y = np.sum(ar_d_mean[:, 0] * ar_d_mean[:, 1])
    x_SR = np.sqrt(np.sum(ar_d_mean[:, 0] * ar_d_mean[:, 0]))
    y_SR = np.sqrt(np.sum(ar_d_mean[:, 1] * ar_d_mean[:, 1]))
    
    corr = round(x_y/(x_SR * y_SR), 2)
    
    return corr

##############################################################################
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
        
    AUC, TPR, FPR = AUC_TPR_FPR(pred, observed_values, 0)
    
    # Calculate the absolute correct
    number_correct = 0; correct_ratio = 0
    if len(selected_cut_DDG) != 0:
        for match in selected_cut_DDG:
            if match[1] * match[2] < 0:
                number_correct += 1
        correct_ratio = number_correct / len(selected_cut_DDG)
    
    return selected_cut_DDG, AUC, TPR, FPR, correct_ratio

'''####################################################################'''

if __name__ == '__main__':  

    preliminary_pred = {}; workable_output = {}
      
    mutation_d = '/home/leo/Documents/Database/Data_Code_Publish/Mutations'
#    mutation_d = '/home/leo/Documents/Database/Pipeline_New/Codes/Results'
    workable_input = Workable_input(mutation_d) # Change the file names inside this function.
    
    os.chdir('/home/leo/Documents/Database/Data_Code_Publish/Structures')
    with open('sequence', 'r') as f:
        sequence = json.load(f)
    with open('combined_ids', 'r') as f:
        combined_ids = json.load(f)
    with open('matched_ids', 'r') as f:
        matched_ids = json.load(f)
    search_para = {}
    search_para['moving'] = True
    search_para['step_size'] = 0.25
    search_para['start_dist'] = 3
    search_para['end_dist'] = 8
    search_para['cut_dist'] = 5.5
    
    for form in ['one']:       
        search_para['form'] = form
        for within_range in [True]:
            search_para['within_range'] = within_range
            
            print('Working on: '+ form +'    '+ str(within_range))
            preliminary_pred[form+'_WithinRange_'+str(within_range)] = {}
            # Container
            container = []

            structure_d = '/home/leo/Documents/Database/Data_Code_Publish/Structures/imgt'    
            workable = Workable_output(workable_input, search_para, combined_ids,matched_ids, sequence, structure_d)
            workable_output[form+'_WithinRange_'+str(within_range)] = workable
            
            working_d = '/home/leo/Documents/Database/Data_Code_Publish/Codes/Results'
            predict_affinity_results = Predict_affinity(workable, working_d, binary= True)
            
            preliminary_pred[form+'_WithinRange_'+str(within_range)]['predict_results_all'] = \
                                            predict_affinity_results
            
            for ran in [[0,0.5], [0, 100],[0.5, 100], [1, 100]]:
                cut_DDG_lower = ran[0]
                cut_DDG_upper = ran[1]
                selected_cut_DDG, AUC, TPR, FPR, correct_ratio = \
                            Analyze_resutls(predict_affinity_results, cut_DDG_lower, cut_DDG_upper)
                concentrations = [Calculate_concentration(TPR, 0.1), Calculate_concentration(TPR, 0.05)]
                container.append([ran, AUC, concentrations[:], len(selected_cut_DDG)])
                
            preliminary_pred[form+'_WithinRange_'+str(within_range)]['range_auc_concentn_len'] =\
                            copy.deepcopy(container)
            
    saving_d = '/home/leo/Documents/Database/Data_Code_Publish/Codes/Results'
#    os.chdir(saving_d)
#    with open('affinity_pre_results', 'w') as f:
#        json.dump(preliminary_pred, f)
#    with open('workable_output', 'w') as f:
#        json.dump(workable_output, f)
#preliminary_pred.keys()
#preliminary_pred['one_WithinRange_True'].keys()
#preliminary_pred['one_WithinRange_True']['range_auc_concentn_len']
#correct_ratio
#preliminary_pred['flanked_WithinRange_False']['predict_results_all']
#predict_affinity_results


'''###################################################################################################'''
'''WE USE THE form multiple_WithinRange_True, and binary model, WHICH IS MORE REASONABLE.'''
          
def Sub_Bootstrap_corr(predictable, iteration):
        # boot
    boot_AUC = []; AUC_CI95 = []; boot_corr=[]; corr_CI95=[]
    for i in range(iteration):
        boot_samples = random.choices(predictable, k = len(predictable))
        # cut the boot sample on basis of differnt DDG
        AUC_unit = []; corr_unit = []
        for ran in [(0, 0.5), (0, 100), (0.5, 100), (1, 100)]:
            DDG_selected_boot = [b for b in boot_samples if abs(b[1])>ran[0] and abs(b[1])<=ran[1]]
            corr = Calculate_correlation(DDG_selected_boot)
            # get the pred and the observed
            pred = [x[2] for x in DDG_selected_boot]
            observed = [x[1] for x in DDG_selected_boot]# use the opposit for the sake of DDG
            
            # Calculate the AUC
            AUC, _, _ = AUC_TPR_FPR(pred, observed, cut = 0)
            # Load
            AUC_unit.append(AUC)
            corr_unit.append(corr)
        # load
        boot_AUC.append(AUC_unit[:])
        boot_corr.append(corr_unit[:])
    # calculate the 95% confidence interval
    for i in range(4):
        # Get the AUC_CI95
        temp_AUC = [x[i] for x in boot_AUC if x[i] != []]# some of the AUC may be []
    
        temp_AUC.sort
        small_ind = math.floor(len(temp_AUC) * 0.025)
        large_ind = math.floor(len(temp_AUC) * 0.975)
    
        AUC_CI95.append([round(temp_AUC[small_ind],2), round(temp_AUC[large_ind],2)])
        
        # Get the corr_CI95
        temp_corr = [x[i] for x in boot_corr]
    
        temp_corr.sort
        small = math.floor(len(temp_corr) * 0.025)
        large = math.floor(len(temp_corr) * 0.975)
        
        corr_CI95.append([round(temp_corr[small],2), round(temp_corr[large],2)])        
    
    return AUC_CI95, corr_CI95
##################################################################################
'''This function has to be carefully tested'''
def Select_other_pred(predictable_mut_ids, original_mut_df, pred_all_df):
    # check the pdbid and the mutation
    selected = []
    for mutid in predictable_mut_ids:
        original_pdb = original_mut_df.iloc[mutid-1]['#PDB']
        original_mut = original_mut_df.iloc[mutid-1]['Mutation'].upper().split(',')
        original_mut.sort()
        for i in range(pred_all_df.shape[0]):
            pred_pdb = pred_all_df.iloc[i]['#pdb'][-4:]
            if pred_pdb == original_pdb:
                pred_mut = pred_all_df.iloc[i]['mutation'].split(',')
                pred_mut.sort()
                if pred_mut == original_mut:
                    selected.append([mutid, pred_all_df.iloc[i]['expt'], pred_all_df.iloc[i]['comp']])
                    break
    return selected
##################################################################################
def Bootstrap_AUC_corr(iteration, WithinRange = True):
    AUC_corr_CI95 = {}
    # Get the pred results of different methods
    os.chdir('/home/leo/Documents/Database/Data_Code_Publish/Codes/Results')
    with open('affinity_pre_results', 'r') as f:        
        affinity_pred = json.load(f)
        
    # Here the key can be changed to get different results under different model
    key = 'multiple_WithinRange_'+str(WithinRange)
    my_pred = affinity_pred[key]['predict_results_all']
    
    # Get rid of the unpredictable
    predictable = []
    for res in my_pred:
        if res[2] != 'Unpredicable' and res[1] != 0:
                predictable.append([res[0], -res[1], res[2]])# Use the opposite of res[1] because of DDG
    # Get the mut ids of the predictable
    predictable_mut_ids = [x[0] for x in predictable]
    # boot
    print('Bootstraping CN')
    AUC_CI95_my, corr_CI95_my = Sub_Bootstrap_corr(predictable, iteration)
    AUC_corr_CI95['CN'] = ('AUC_CI95:', AUC_CI95_my, 'corr_CI95', corr_CI95_my)
        # Boot other results
    os.chdir('/home/leo/Documents/Database/Data_Code_Publish/Mutations')
    original_mutations = pd.read_excel('Mutation.xlsx')
    # Read the original mutation 
    for method in ['bASA', 'dDfire','dfire', 'discovery_studio', 'rosetta', 'statium', 'foldX']:
        print('Bootstraping '+ method)        
        scores = pd.read_csv(method+'.dat', sep = ';')    
        selected = Select_other_pred(predictable_mut_ids, original_mutations, scores)
        AUC_CI95_other, corr_CI95_other = Sub_Bootstrap_corr(selected , iteration)
        
        AUC_corr_CI95[method] = ('AUC_CI95:', AUC_CI95_other, 'corr_CI95', corr_CI95_other)
  
    return AUC_corr_CI95

#########################################################################
#CI95 = {}
#for WithinRange in [True, False]:    
#    AUC_corr_CI95 = Bootstrap_AUC_corr(10000)
#    CI95['multiple_WithinRange_'+str(WithinRange)] = AUC_corr_CI95
#   
#os.chdir('/home/leo/Documents/Database/Data_Code_Publish/Codes/Results')
#
#with open('AUC_corr_CI95', 'r') as f:  
#    CI95_read = json.load(f)
#CI95_read  


'''Change the file name in the Workable function can change the source of mutations'''

# CHECK IF THE WORKABLE IS CONSISTANT WITH THE GIVEN
#mutation_d = '/home/leo/Documents/Database/Data_Code_Publish/Mutations'
#
#workable_input = Workable_input(mutation_d) 
#
#workable_input
#
#search_para = {}
#search_para['moving'] = True
#search_para['step_size'] = 0.25
#search_para['start_dist'] = 3
#search_para['end_dist'] = 8
#search_para['cut_dist'] = 6
#
#search_para['within_range'] = True
#search_para['form'] = 'one'
#
#
#structure_d = '/home/leo/Documents/Database/Data_Code_Publish/Structures/imgt' 
#workable_output = \
#        Workable_output(workable_input, search_para, combined_ids, sequence, structure_d)
#
#
#mut_id = 90
#
#len(workable_output)
#len(workable_input)
#
#lower  = 532
#upper = 537
#
#workable_output[lower:upper]
#workable_input[lower:upper]
#
#combined_ids['3be1']



#
#
#os.chdir('/home/leo/Documents/Database/Pipeline_New/Codes/Results')
#with open('one_to_one_results', 'r') as f:
#    one_to_one = json.load(f)
#
#
##one_to_one.keys()
##one_to_one['affinity_results']['results_dict'].keys()
##one_to_one['affinity_results']['results_dict']['single_True_mode_single_moving_True_binary_True'].keys()
#res = one_to_one['affinity_results']['results_dict']['single_True_mode_single_moving_True_binary_True']['0']
#
#old_workable = []
#for complicated in res:
#    mutation = complicated[1:len(complicated)-1]
#
#    mut_id = complicated[0][3]
#    DDG = complicated[0][2][2]
#
#    mu_set = []
#    for mmu in mutation:
#        if mmu != []:
#            mu_set.append([mmu[0][0], mmu[1][0], 1, DDG, mut_id])
#    old_workable.append(mu_set[:])
#
#
#working_d = '/home/leo/Documents/Database/Data_Code_Publish/Codes/Results'
#predict_affinity_results = Predict_affinity(old_workable, working_d, binary = True)
#container = []
#for ran in [[0,0.5], [0, 100],[0.5, 100], [1, 100]]:
#    cut_DDG_lower = ran[0]
#    cut_DDG_upper = ran[1]
#    selected_cut_DDG, AUC, TPR, FPR, correct_ratio = \
#                Analyze_resutls(predict_affinity_results, cut_DDG_lower, cut_DDG_upper)
#    concentrations = [Calculate_concentration(TPR, 0.1), Calculate_concentration(TPR, 0.05)]
#    container.append([ran, AUC, concentrations[:], len(selected_cut_DDG)])
#    
#container



