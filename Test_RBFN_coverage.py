'''CODE FOR TESTING'''
import numpy as np
import os
import json
import math
import copy
from matplotlib import pyplot as plt

###########################################################
os.chdir('/home/leo/Documents/Database/Data_Code_Publish/Codes')
from AlignmentConstraint import To_seq
from Train_Loss import Train
from RBFN_coverage import Distance_matrix, Design_matrix, Observed_values, Calculate_AUC
#############################################################################
# Define  distances
from Bio import Align
from Bio.SubsMat.MatrixInfo import blosum62
aligner = Align.PairwiseAligner()
aligner.substitution_matrix = blosum62
aligner.mode = 'global'

#############################################################################
def Test_discrimination(train_results, positive_d, negative_pre):
    test_discrimination = {}                
    # Extract the observed values
    for binary in ['binary', 'numerical']:        
        # Take out the gap_score and the extened gap score
        key = 'cross_'+ binary + '_Gaussian_'
        best_gap_ext = train_results[key]['best_gap_ext']
        # Pass to the aligner
        aligner.gap_score = best_gap_ext[0]
        aligner.extend_gap_score = best_gap_ext[1]
        
        test_discrimination[key] = {}
        
        for i in range(1, 5):
            for j in range(1,5):
                k = str(i)+'_'+str(j)+'_0_0_1_2_1perchain'
                # Take out the cneters and the coefficients
                centers = train_results[key][k]['centers']
                coeff = train_results[key][k]['coefficients']
                
                # Prepare the testing data
                p_name = 'testing_' + k
                n_name = p_name + '_negative'
                
                # Open the positive testing set
                os.chdir(positive_d)
                with open(p_name, 'r') as f:
                    positive_test = json.load(f)
                    
                test_discrimination[key][k]={}
                test_discrimination[key][k]['AUC_list'] = []
                test_discrimination[key][k]['TPR_list'] = []
                test_discrimination[key][k]['FPR_list'] = []
                    
                # Different negative set
                for n in range(10):
                    negative_d = negative_pre + str(n)
                    os.chdir(negative_d)
                    with open(n_name, 'r') as f:
                        negative_test = json.load(f)
                        
                    #Prepare the testing set
                    testing_set = copy.deepcopy(positive_test)
                    testing_set.extend(negative_test)
                    
                    # Calculate the observed values
                    observed_test = Observed_values(testing_set, binary)
                    
                    # Calculate the distance matrix
                    distance_matrix_test = Distance_matrix(testing_set, centers, square=False)
                    
                    # Calculate the design matrix
                    radius_coeff = np.ones((distance_matrix_test.shape[1], 1))
                    design_matrix_test = Design_matrix(distance_matrix_test,radius_coeff, basis_function = 'Gaussian')
                    design_matrix_test = np.hstack((design_matrix_test, np.ones((len(testing_set), 1))))
                    
                    # Make prediction
                    coeff = np.reshape(coeff, (-1,1))
                    pred = design_matrix_test.dot(coeff)
                    
                    # Calculate the AUC, FPR, TRP
                    AUC, TPR, FPR = Calculate_AUC(pred, observed_test)
                    print(AUC)
                    # Load the results
                    test_discrimination[key][k]['AUC_list'].append(AUC)
                    test_discrimination[key][k]['TPR_list'].append(TPR)
                    test_discrimination[key][k]['FPR_list'].append(FPR)
                    
    return test_discrimination
'''#################################################################################3'''
def Test_reverse(train_results, positive_d):

    test_reverse = {}                
    # Extract the observed values
    for binary in ['binary', 'numerical']:        
        # Take out the gap_score and the extened gap score
        key = 'cross_'+ binary + '_Gaussian_'
        best_gap_ext = train_results[key]['best_gap_ext']
        # Pass to the aligner
        aligner.gap_score = best_gap_ext[0]
        aligner.extend_gap_score = best_gap_ext[1]
        
        test_reverse[key] = {}
        
        for i in range(1, 5):
            for j in range(1,5):
                k = str(i)+'_'+str(j)+'_0_0_1_2_1perchain'
                # Take out the cneters and the coefficients
                centers = train_results[key][k]['centers']
                coeff = train_results[key][k]['coefficients']
                
                # Prepare the testing data
                p_name = 'testing_' + k
                n_name = 'testing_'+str(j)+'_'+str(i)+'_0_0_1_2_1perchain'
                
                # Open the positive testing set
                os.chdir(positive_d)
                with open(p_name, 'r') as f:
                    positive_test = json.load(f)
                with open(n_name, 'r') as f:
                    negative_test_reverse = json.load(f)
                
                # Reverse the order of the negative_test and set the oberved to 0
                negative_test = []
                for core in negative_test_reverse:
                    negative_test.append([core[1], core[0], 0, 0])
                    
                              
                    
                #Prepare the testing set
                testing_set = copy.deepcopy(positive_test)
                testing_set.extend(negative_test)
                
                # Calculate the observed values
                observed_test = Observed_values(testing_set, binary)
                
                # Calculate the distance matrix
                distance_matrix_test = Distance_matrix(testing_set, centers, square=False)
                
                # Calculate the design matrix
                radius_coeff = np.ones((distance_matrix_test.shape[1], 1))
                design_matrix_test = Design_matrix(distance_matrix_test,radius_coeff, basis_function = 'Gaussian')
                design_matrix_test = np.hstack((design_matrix_test, np.ones((len(testing_set), 1))))
                
                # Make prediction
                coeff = np.reshape(coeff, (-1,1))
                pred = design_matrix_test.dot(coeff)
                
                # Calculate the AUC, FPR, TRP
                AUC, TPR, FPR = Calculate_AUC(pred, observed_test)
                print(AUC)
                # Load the results
                test_reverse[key][k]={} 
                test_reverse[key][k]['AUC']=AUC
                test_reverse[key][k]['TPR'] = TPR
                test_reverse[key][k]['FPR'] = FPR
    
    return test_reverse
    

if __name__ == '__main__':
    os.chdir('/home/leo/Documents/Database/Data_Code_Publish/Codes/Results')
    with open('train_results', 'r') as f:
        train_results = json.load(f)
    positive_d = '/home/leo/Documents/Database/Data_Code_Publish/Cores/Positive_cores'
    negative_pre = '/home/leo/Documents/Database/Data_Code_Publish/Cores/Negative_cores/Sample_'
    
    test_discrimination = Test_discrimination(train_results, positive_d, negative_pre)
    test_reverse = Test_reverse(train_results, positive_d)

    # Save the results
    saving_d = '/home/leo/Documents/Database/Data_Code_Publish/Codes/Results'
    os.chdir(saving_d)
    with open('test_discrimination', 'w') as f:
        json.dump(test_discrimination, f)        
    with open('test_reverse', 'w') as f:
        json.dump(test_reverse, f)























