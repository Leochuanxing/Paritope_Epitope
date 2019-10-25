'''THIS MODULE IS AN INTEGRATED MODEL'''  

import json
import os
import numpy as np

from Bio import Align
from Bio.SubsMat.MatrixInfo import blosum62
aligner = Align.PairwiseAligner()
aligner.substitution_matrix = blosum62
aligner.mode = 'global'


os.chdir('/home/leo/Documents/Database/Data_Code_Publish/Codes')
from RBFN_coverage import Distance_matrix, Design_matrix
from Some_basic_fuctions import AUC_TPR_FPR

'''#####################################################################'''        
'''
Input: matched pair
Output: a score
'''         
 
       
def One_integrated(matched_pair, binary, working_d):
    Ab_aa = matched_pair[0]
    Ag_aa = matched_pair[1]
    cn = matched_pair[2]
    # Generate combinations
    one_one_combination=[]
    for ab in Ab_aa:
        for ag in Ag_aa:
            one_one_combination.append([[ab], [ag]])
            
    os.chdir(working_d)
    with open('train_results', 'r') as f:
        train_results = json.load(f)
        
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
    
    key = '1_1_0_0_1_2_1perchain'
    # Get the centers and coefficients for related match type
    centers = model_all[key]['centers']
    coefficients = model_all[key]['coefficients']
    
    # Calculate the design matrix

    distance_matrix = Distance_matrix(one_one_combination, centers, square=False)
    radius_coeff = np.ones((len(centers), 1))
    design_matrix = Design_matrix(distance_matrix, radius_coeff)
    # Attach the constant column
    design_matrix = np.hstack((design_matrix, np.ones((len(one_one_combination), 1))))
    
    # Predict
    coefficients = np.reshape(coefficients, (-1, 1))
    pred = np.average(design_matrix.dot(coefficients))
    
    # pred_observed
    pred_observed = [pred, cn]
    
    return pred_observed
'''########################################################################3'''       
def Integrated(match_type_set, binary, working_d):
    pred =[]; observed = []
    for matched_pair in match_type_set:
        pred_observed = One_integrated(matched_pair, binary, working_d)
        pred.append(pred_observed[0])
        observed.append(pred_observed[1])
    
    # Calculate auc
    AUC, TPR, FPR = AUC_TPR_FPR(pred, observed, 0.1)
        
    return AUC, TPR, FPR
'''######################################################################'''
def Integrated_discrimination():
    pass
'''######################################################################'''
if __name__ == '__main__':
    working_d = '/home/leo/Documents/Database/Data_Code_Publish/Codes/Results'
    integrated_results = {}; integrated_reverse = {}
    for i in range(1,4):
        for j in range(1, 4):
            positive_name = 'testing_'+str(i)+'_'+str(j)+'_0_0_1_2_1perchain'
            print('Working on  ', positive_name)
            integrated_results[positive_name] = {}; integrated_reverse[positive_name]={}
            os.chdir('/home/leo/Documents/Database/Data_Code_Publish/Cores/Positive_cores')
            with open(positive_name, 'r') as f:
                positive = json.load(f)
            auc_list = []; tpr_list = []; fpr_list = []
            for k in range(10):
                os.chdir('/home/leo/Documents/Database/Data_Code_Publish/Cores/Negative_cores/Sample_'+str(k))
                negative_name = positive_name + '_negative'
                with open(negative_name, 'r') as f:
                    negative = json.load(f)
                
                testing = []
                testing.extend(positive)
                testing.extend(negative)
                
                AUC, TPR, FPR = Integrated(testing, binary=True, working_d=working_d)
                auc_list.append(AUC)
                tpr_list.append(TPR)
                fpr_list.append(FPR)
                
            integrated_results[positive_name]['AUC_list'] = auc_list
            integrated_results[positive_name]['TPR_list'] = tpr_list
            integrated_results[positive_name]['FPR_list'] = fpr_list
            
            # Calculate the reverse
            reverse_name = 'testing_'+str(j)+'_'+str(i)+'_0_0_1_2_1perchain'
            os.chdir('/home/leo/Documents/Database/Data_Code_Publish/Cores/Positive_cores')
            with open(reverse_name, 'r') as f:
                reverse = json.load(f)
            reverse_negative = []
            for match in reverse:
                reverse_negative.append([match[1], match[0], 0])
            
            testing_reverse = []
            testing_reverse.extend(positive)
            testing_reverse.extend(reverse_negative)
            AUC, TPR, FPR = Integrated(testing_reverse, binary=True, working_d=working_d)
            
            integrated_reverse[positive_name]['AUC'] = AUC
            integrated_reverse[positive_name]['TPR'] = TPR
            integrated_reverse[positive_name]['FPR'] = FPR
            
        
        
            
              
    os.chdir('/home/leo/Documents/Database/Data_Code_Publish/Codes/Results') 
    with open('integrated_discrimination', 'w') as f:
        json.dump(integrated_results, f)    
    with open('integrated_reverse', 'w') as f:
        json.dump(integrated_reverse, f)
        
    
                

    
    
    
    
    
    
    
    