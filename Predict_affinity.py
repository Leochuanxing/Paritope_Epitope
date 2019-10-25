'''PREDICT MUTATION'''
import sys
import json
import math
import os
import random
import numpy as np
import pandas as pd

'''#############################################################'''
from Bio import Align
from Bio.SubsMat.MatrixInfo import blosum62
aligner = Align.PairwiseAligner()
aligner.substitution_matrix = blosum62
aligner.mode = 'global'
'''###################################################################'''
os.chdir('/home/leo/Documents/Database/Data_Code_Publish/Codes')
from RBFN_coverage import Distance_matrix, Design_matrix
from RBFN_coverage import   NumpyEncoder
from Some_basic_fuctions import AUC_TPR_FPR
'''######################################################################'''
'''

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

'''
Input:
    cut_DDG_lower, cut_DDG_upper: cut values, only to consider the mutations with 
            the DDG in 
                        (cut_DDG_lower, cut_DDG_upper]
            open in the left and closed in the right.
'''
def Analyze_resutls(predict_affinity_results, cut_range):
    # Select the prediction by cut_DDG
    AUC_list = []; TPR_list = []; FPR_list = []
    for ran in cut_range:
        cut_DDG_lower = ran[0]
        cut_DDG_upper = ran[1]
        selected_cut_DDG = []
        for res in predict_affinity_results:
            if res[2] != 'Unpredicable' and res[1] != 0:
                if abs(res[1]) > cut_DDG_lower and abs(res[1]) <= cut_DDG_upper:
                    selected_cut_DDG.append(res)
        # Calculate AUC, FPR, TPR
        observed_values = []; pred = []
        for match in selected_cut_DDG:
            # We have to add a nagetive sign to the observed to make sure the 
            # observed_values and the pred have monotonicity
            observed_values.append(match[1])
            pred.append(match[2])
            
        AUC, TPR, FPR = AUC_TPR_FPR(pred, observed_values, 0)
        
        AUC_list.append(AUC)
        TPR_list.append(TPR)
        FPR_list.append(FPR)
    
    return  AUC_list, TPR_list, FPR_list

'''####################################################################'''

          
def Sub_Bootstrap_corr(predictable, cut_range, iteration):
        # boot
    boot_AUC = []; AUC_CI95 = []
    for i in range(iteration):
        boot_samples = random.choices(predictable, k = len(predictable))
        # cut the boot sample on basis of differnt DDG
        AUC_unit = []
        for ran in cut_range:
            DDG_selected_boot = [b for b in boot_samples if abs(b[1])>ran[0] and abs(b[1])<=ran[1]]

            # get the pred and the observed
            pred = [x[2] for x in DDG_selected_boot]
            observed = [x[1] for x in DDG_selected_boot]
            
            # Calculate the AUC
            AUC, _, _ = AUC_TPR_FPR(pred, observed, cut = 0)
            # Load
            AUC_unit.append(AUC)

        # load
        boot_AUC.append(AUC_unit[:])

    # calculate the 95% confidence interval
    for i in range(len(cut_range)):
        # Get the AUC_CI95
        temp_AUC = [x[i] for x in boot_AUC if x[i] != []]# some of the AUC may be []
    
        temp_AUC.sort(reverse = False)
        small_ind = math.floor(len(temp_AUC) * 0.025)
        large_ind = math.floor(len(temp_AUC) * 0.975)
    
        AUC_CI95.append([round(temp_AUC[small_ind],2), round(temp_AUC[large_ind],2)])
             
    
    return AUC_CI95
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


#########################################################################
if __name__ == '__main__':
    pred_workable_formated = {}
    
    cut_range = [(0, 0.5), (0, 100), (0.5, 100), (1, 100)]
    iteration = 10000 # the iteration number for boot straping
    
    pred_workable_formated['cut_range'] = cut_range
    pred_workable_formated['iteration'] = iteration
    # Load the data
    working_d = '/home/leo/Documents/Database/Data_Code_Publish/Codes/Results'
    os.chdir(working_d)
    with open('workable_formated', 'r') as f:
        workable_formated = json.load(f)
    # The results under different mode
    for key, workable in workable_formated.items():
#    key = 'one_WithinRange_True'
#    workable = workable_formated[key]
        print('Working one mode    ', key)
        pred_workable_formated[key]={}
        # Get the pred results
        pred_results = Predict_affinity(workable, working_d, binary = True)
        pred_workable_formated[key]['pred_results'] = pred_results
        
        # Extract the predictable
        predictable = []
        for res in pred_results:
            if res[2] != 'Unpredicable' and res[1] != 0:
                    predictable.append([res[0], -res[1], res[2]])# Use the opposite of res[1] because of DDG
        
        # Calculate AUC, FPR, TPR CI95,corr, concentration and draw AUROC
        # CN method
        AUC_list, TPR_list, FPR_list = Analyze_resutls(predictable, cut_range)
        pred_workable_formated[key]['CN'] = {}
        pred_workable_formated[key]['CN']['AUC_list'] = AUC_list[:]
        pred_workable_formated[key]['CN']['TPR_list'] = TPR_list[:]
        pred_workable_formated[key]['CN']['FPR_list'] = FPR_list[:]
        print('Bootstraping CN')
        AUC_CI95_my = Sub_Bootstrap_corr(predictable, cut_range, iteration)
        pred_workable_formated[key]['CN']['BootStrap'] = ('AUC_CI95:', AUC_CI95_my)        
        
        # Other methods
        # Get the mut ids of the predictable
        predictable_mut_ids = [x[0] for x in predictable]
        os.chdir('/home/leo/Documents/Database/Data_Code_Publish/Mutations')
        original_mutations = pd.read_excel('Mutation.xlsx')
        # Read the original mutation 
        for method in ['bASA', 'dDfire','dfire', 'discovery_studio', 'rosetta', 'statium', 'foldX']:       
            scores = pd.read_csv(method+'.dat', sep = ';')    
            selected = Select_other_pred(predictable_mut_ids, original_mutations, scores)
            AUC_list, TPR_list, FPR_list = Analyze_resutls(selected, cut_range)
            pred_workable_formated[key][method] = {}
            pred_workable_formated[key][method]['AUC_list'] = AUC_list[:]
            pred_workable_formated[key][method]['TPR_list'] = TPR_list[:]
            pred_workable_formated[key][method]['FPR_list'] = FPR_list[:]
            print('Bootstraping '+ method)
            AUC_CI95_other= Sub_Bootstrap_corr(selected,cut_range, iteration)
            pred_workable_formated[key][method]['BootStrap'] = ('AUC_CI95:', AUC_CI95_other)
    # Save the results
#    os.chdir(working_d)
#    with open('pred_workable_skempi', 'w') as f:
#        json.dump(pred_workable_formated, f, cls=NumpyEncoder)
        
#for me in ['CN', 'bASA', 'dDfire','dfire', 'discovery_studio', 'rosetta', 'statium', 'foldX']: 
#    print(me)     

#    print(pred_workable_formated['flanked_WithinRange_True'][me]['AUC_list'])
###    print(pred_workable_formated['flanked_WithinRange_True'][me]['AUC_list'])
#    print([len(tpr) for tpr in pred_workable_formated['flanked_WithinRange_True'][me]['TPR_list']])
#    print(pred_workable_formated['flanked_WithinRange_True'][me]['BootStrap'] )
#    print('n/')


