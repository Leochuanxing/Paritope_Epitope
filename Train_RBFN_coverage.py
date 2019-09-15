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
from RBFN_coverage import Distance_matrix, Design_matrix, Observed_values
from RBFN_coverage import  Coverage_reduce_centers, NumpyEncoder

#############################################################
# Define  distances
from Bio import Align
from Bio.SubsMat.MatrixInfo import blosum62
aligner = Align.PairwiseAligner()
aligner.substitution_matrix = blosum62
aligner.mode = 'global'
###########################################################
def Train_RBFN_coverage(best_parameter,positive_d, negative_d):
    # The training results will be stored inn train_results
    train_results = {}                
    # Extract the observed values
    for binary in ['binary', 'numerical']:        
        # Take out the gap_score and the extened gap score
        key = 'cross_'+ binary + '_Gaussian_'
        best_gap_ext = best_parameter[key]['best_gap_ext']
        # Pass the parameter to the aligner
        aligner.gap_score = best_gap_ext[0]
        aligner.extend_gap_score = best_gap_ext[1]
        # Pass the parameter to the train_results
        train_results[key] = {}
        train_results[key]['best_gap_ext'] = best_gap_ext
        
        for i in range(1, 5):
            for j in range(1,5):
                
                # Open the positive cores and negative cores
                p_name = 'training_'+str(i)+'_'+str(j)+'_0_0_1_2_1perchain'
                os.chdir(positive_d)
                with open(p_name, 'r') as f:
                    positive_training_set = json.load(f)
                n_name = p_name+'_negative'
                os.chdir(negative_d)
                with open(n_name, 'r') as f:
                    negative_training_set = json.load(f)
                
                # Take out the best reg and percentage
                k = str(i)+'_'+str(j)+'_0_0_1_2_1perchain'
                reg = best_parameter[key][k][0]
                percentage = best_parameter[key][k][1]
                
                print(p_name)
                print('Calculating the distance matrix and the design matrix')
                # Calculate the distance matrices and the design matrices
                training_set = copy.deepcopy(positive_training_set)
                training_set.extend(negative_training_set)
                # Calculate the distance matrix
                all_training_distance_matrix = Distance_matrix(training_set, training_set, square = True)
                # Calculate the design matrix
                radius_coeff = np.ones((all_training_distance_matrix.shape[1], 1))
                all_training_design_matrix = Design_matrix(all_training_distance_matrix,radius_coeff, basis_function='Gaussian')
                
                # Select the centers according to the percentage
                # Get the beginning center indices
                centers_init = []; non_redundant = []
                for indx in range(len(training_set)):
                    if [training_set[indx][0],training_set[indx][1]] not in non_redundant:
                        non_redundant.append([training_set[indx][0],training_set[indx][1]])
                        centers_init.append(indx)
                # Extract the distance matrix for selecting the centers
                distance_matrix_center = all_training_distance_matrix[centers_init,:][:, centers_init]
                # Calculate the number of centers
                n_centers = math.floor(len(centers_init) * percentage)
                # Select the center indices corresponding to training_set
                centers_ind_list = Coverage_reduce_centers(distance_matrix_center,centers_init,[n_centers])
                centers_ind = centers_ind_list[0]# Those indices can are corresponding to the training_set
                
                # Extract the design matrix
                design_matrix = all_training_design_matrix[:, centers_ind]
                design_matrix = np.hstack((design_matrix, np.ones((len(training_set), 1))))
                
                # Calculate the observed values
                observed_values = Observed_values(training_set, binary)
                
                # Load the values to train_para
                train_para = {}
                train_para['observed'] = observed_values.reshape(-1, 1)
                train_para['design_matrix'] = design_matrix
                train_para['reg'] = reg
                train_para['method'] = 'BFGS'
#                train_para['step_size']
#                train_para['n_iterations']
                
                if binary:
                    train_para['loss_type'] = 'Sigmoid'
                else:
                    train_para['loss_type'] = 'SumSquares'
                
                coeff = np.zeros((design_matrix.shape[1], 1))
                train_para['coefficients'] = coeff.reshape(-1, 1)
                
                # Train the model
                print('Training')
                termination = 1E-4 * n_centers
                train_para, loss = Train(train_para, rho=0.8, c = 1e-4, termination = termination)
                
                # Store and save the training results
                train_results[key][k] = {}
                train_results[key][k]['coefficients'] = copy.deepcopy(train_para['coefficients'])
                train_results[key][k]['centers'] = [training_set[x] for x in centers_ind]
                
    return train_results

'''###############################################################################'''

if __name__ == '__main__':
    
    positive_d = '/home/leo/Documents/Database/Data_Code_Publish/Cores/Positive_cores'
    negative_d = '/home/leo/Documents/Database/Data_Code_Publish/Cores/Negative_cores/Sample_0'
    saving_d = '/home/leo/Documents/Database/Data_Code_Publish/Codes/Results'
    
    best_parameter_d = '/home/leo/Documents/Database/Data_Code_Publish/Codes/Cross_validation'
    os.chdir(best_parameter_d)
    with open('best_parameter', 'r') as f:
        best_parameter = json.load(f)
    
    train_results = Train_RBFN_coverage(best_parameter,positive_d, negative_d)
    
    # Save the train results
    os.chdir(saving_d)
    with open('train_results', 'w') as f:
        json.dump(train_results, f, cls=NumpyEncoder)


#os.chdir(saving_d)
#with open('train_results', 'r') as f:
#    train_results = json.load(f)
#train_results.keys()
#train_results['cross_binary_Gaussian_'].keys()
#train_results['cross_binary_Gaussian_']['1_1_0_0_1_2_1perchain'].keys()
#type(train_results['cross_binary_Gaussian_']['1_1_0_0_1_2_1perchain']['coefficients'])
#type(train_results['cross_binary_Gaussian_']['1_1_0_0_1_2_1perchain']['centers'])
#len(train_results['cross_binary_Gaussian_']['1_1_0_0_1_2_1perchain']['coefficients'])
#len(train_results['cross_binary_Gaussian_']['1_1_0_0_1_2_1perchain']['centers'])
##    core = json.load(f)
#train_results['cross_binary_Gaussian_']['1_1_0_0_1_2_1perchain']['coefficients'][-3:]
#train_results['cross_binary_Gaussian_']['1_1_0_0_1_2_1perchain']['centers'][:4]
#best_parameter





































