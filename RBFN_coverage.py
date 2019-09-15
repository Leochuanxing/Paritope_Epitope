###############################################################
# Import the modules

#import random
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

#############################################################
# Define  distances
from Bio import Align
from Bio.SubsMat.MatrixInfo import blosum62
aligner = Align.PairwiseAligner()
aligner.substitution_matrix = blosum62
aligner.mode = 'global'
################################################################


'''
Multiplication_distance:
    This function is to calculate the multiplication distance for a two given matched paires
Input:
    parepi_1:
        a matched paire, in the form of [['ALA', 'SER'], ['CYS', 'GLU'],...]
        ['ALA', 'SER'] are from the antibody, and ['CYS', 'GLU'] are from the antigen.
        parepi_1 can be of any length as long as the first two elements are as required.
    parepi_2:
        in the same form as parepi_1
Output:
    distance:
        a float, gives the distance between parepi_1 and parepi_2
'''
def Multiplication_distance(parepi_1, parepi_2):
    Ab_seq1 = To_seq(parepi_1[0])
    Ab_seq2 = To_seq(parepi_2[0])
    Ag_seq1 = To_seq(parepi_1[1])
    Ag_seq2 = To_seq(parepi_2[1])
    l_Ab = len(Ab_seq1)
    l_Ag = len(Ag_seq1)
    distance = (1 - (4*l_Ab + aligner.score(Ab_seq1, Ab_seq2))/(15*l_Ab)) *  (1 - (4*l_Ag + aligner.score(Ag_seq1, Ag_seq2))/(15*l_Ag))
    return distance  

################################################################################
'''
Distance_matrix:
    this function is to calculate the distance matrix for the give row-set and col_set. The function 
    used to calculate the distance between the element of the row_set and the element of the col_set 
    is 'Multiplication_distance'.
Input:
    row_set:
        a set of matched paires which decide the rows of the matrix. the form of the 
        row set is given as
       [[['ALA', 'SER'], ['CYS', 'GLU'],....],...]
       Its elements are in the form of [['ALA', 'SER'], ['CYS', 'GLU'],....], 
    col_set:
        a set of matched paires which decide the columns of the distance matrix. It is 
        in the same form as the row_set.
    square:
        it takes the value of either True or False
        When it takes the value of True, it means the row_set and the col_set are the same.
        When it takes the value of False, it means the row_set and the col_set are different.
Output:
    distance_matrix
    
'''
def Distance_matrix(row_set, col_set, square):
    distance_matrix = np.zeros((len(row_set), len(col_set)))
    if square:
        for i in range(len(row_set)):
            for j in range(i, len(row_set)):
                distance_matrix[i, j] = Multiplication_distance(row_set[i], row_set[j])
                distance_matrix[j,i] = distance_matrix[i,j]
    elif not square:
        for i in range(len(row_set)):
            for j in range(len(col_set)):
                distance_matrix[i,j]=Multiplication_distance(row_set[i], col_set[j])
    return distance_matrix
################################################################################
'''
Some basis functions
'''
def Gaussian(distance, radius):
    return math.exp(-(distance**2)*(radius**2))
    
def Mrakov(distance, radius):
    return math.exp(-distance*(radius**2))
'''
beta is a positive number
'''
def Inverse_Multi_Quadric(distance,c, beta):
    return (distance**2 + c**2)**(-beta)
################################################################################
'''
Design_matrix:
    a function to calculate the design matrix according to the distance matrix.
    In our RBFN, we only use 'Markov' as basis function
Input: radius, an array with dimension (n_col, 1), where the n_col is the number of
     columns of the distance_matrix
'''
def Design_matrix(distance_matrix, radius_coeff, basis_function = 'Gaussian'):
    nrow = np.shape(distance_matrix)[0]
    ncol = np.shape(distance_matrix)[1]
    
    design_matrix = np.zeros_like(distance_matrix)
    for i in range(nrow):
        for j in range(ncol):
            if basis_function == 'Gaussian':
                design_matrix[i, j] = Gaussian(distance_matrix[i, j], radius_coeff[j,0])
            elif basis_function == 'Markov':
                design_matrix[i, j] = Mrakov(distance_matrix[i, j], radius_coeff[j,0])
            elif basis_function == 'Inverse_Multi_Quadric':
                design_matrix[i, j] = Inverse_Multi_Quadric(distance_matrix[i, j], c=1, beta=2)  
                
    # Add a constant term to the design matrix
#    design_matrix = np.hstack((design_matrix, np.ones((nrow,1))))
           
    return design_matrix
###################################################################################

'''
Pruning:
    to return the center indices 
Input: 
      training_training_distance_matrix: 
          the distance matrix of all the training samples.
      n_centers_list:
          a list gives the number of centers to be kept
          
***********************The training_training_distance_matrix*****************************
It can be the non_redundant_training_set, and the returned centers_list corresponding to this training
set, which can be connected to the training set through the centers generate by the function 
'Remove_duplicates'. This training_training_distance_matrix should be a square matrix.

Output: 
    centers_list:
        A list of center indices with the length of each list corresponding to n_centers_list
'''

def Coverage_reduce_centers(training_training_distance_matrix, training_set,n_centers_list):
    # Change the main diagonal of the training_training_distance_matrix into 0s
    copy_matrix = copy.deepcopy(training_training_distance_matrix)
    for i in range(np.shape(copy_matrix)[0]):
        copy_matrix[i,i] = 0
        
    from scipy.cluster.hierarchy import linkage, cut_tree
    from scipy.spatial.distance import squareform
    hcluster_CDRH_unordered = linkage(squareform(copy_matrix), method = 'complete')
    
    # Cut the tree 
    cuts = []
    for n in n_centers_list:
        cut_array = cut_tree(hcluster_CDRH_unordered, n_clusters=n)
        cut_list = cut_array.flatten().tolist()
        cuts.append(copy.deepcopy(cut_list))
        
    # Change the cuts into indices, in each cluster, we choose the smallest index as the representive
    centers_list=[]
    for cut in cuts:
        centers_indx = []
        cut_set = set(cut)
        for ind in range(len(cut)):
            if cut[ind] in cut_set:
                centers_indx.append(ind)
                cut_set.remove(cut[ind])
        
        # Find the corresponding elements in the training_set
        centers = []       
        for indx in centers_indx:
            centers.append(training_set[indx])
        centers_list.append(copy.deepcopy(centers))
        
    return centers_list

#####################################################################################
'''
'''
'''
Generate_cross_testing_training:
    A function, generates the testing and training positions for the given number of 
    cross validation. The return values can be used to extract the design matrix from
    the big design matrix
Input:
    training_set:
        a list contains all the training set
    cross_number:
        an integer, gives the number of crosses
Output:
    cross_testing_set:
        a list contains 'cross_number' of lists of positions of the testing_set
    cross_training_set:
        a list contains 'cross_number' of lists of positions of the training_set
        
    The above two list should be corresponding with each other, that is they are 
    ordered correspondingly.
    
*************** A simple example*************************
training_set = [1,2,3,4,5]
cross_number = 2
cross_training_set = [[1,2], [3,4,5]]
cross_testing_set = [[3,4,5], [1,2]]
'''
def Generate_cross_testing_training(training_set, cross_number):
    cross_training_set = []
    cross_testing_set = []
    n_total = len(training_set)
    size = math.floor(n_total/cross_number)
    # Cut into pieces
    for i in range(cross_number-1):
        lower=i*size; upper = (i+1) * size
        cross_testing_set.append(list(range(lower, upper)))
    # Deal with the last block
    lower = size * (cross_number-1)
    cross_testing_set.append(list(range(lower, n_total)))
    # Load the cross_training_set
    for one_test_set in cross_testing_set:
        one_train_set = []
        for j in range(n_total):
            if j not in one_test_set:
                one_train_set.append(j)
        cross_training_set.append(one_train_set)   
            
    return cross_testing_set, cross_training_set

#a = list(range(10))
#cross_test, cross_train = Generate_cross_testing_training(a, 3)
#cross_test
#cross_train
######################################################################################
'''
Raised_cross_indices:
    This function is to generate the cross_training_indices and the cross_testing_indices, so that they 
    can be used directly in the set of 'positive_training_set.extend(negative_training_set)'
Input:
    Positive_training_set:
        a list of positive samples
    Negative_training_set:
        a list of negative samples
    cross_number
Output:
    cross_train_indices: A list of indices correspondint to the design_matrix, observed_values, and the training_set
    cross_test_indices: the same as above
'''
def Raised_cross_indices(positive_training_set, negative_training_set, cross_number):
    n_positive = len(positive_training_set)
    # Generate the cross sets
    positive_cross_test, positive_cross_train = Generate_cross_testing_training(positive_training_set, cross_number)
    negative_cross_test, negative_cross_train = Generate_cross_testing_training(negative_training_set, cross_number)
    
    # Raise the indices of the negative_cross_test and the negative_cross_train by 'n_positive']
    for i in range(cross_number):
        raised_negative_train = (np.array(negative_cross_train[i]) + n_positive).flatten().tolist()
        raised_negative_test = (np.array(negative_cross_test[i]) + n_positive).flatten().tolist()
        positive_cross_train[i].extend(raised_negative_train)
        positive_cross_test[i].extend(raised_negative_test)
    
    cross_train_indices = positive_cross_train
    cross_test_indices = positive_cross_test
    
    return cross_train_indices, cross_test_indices
#a = [1, 2, 3, 4, 5, 6, 7]
#b = [5, 4, 6, 4, 5, 5, 5]
#cross_train_indices, cross_test_indices = Raised_cross_indices(a, b, 3)
#cross_train_indices
#cross_test_indices
###########################################################################################
'''
The purpose of this block is to do cross validation, select centers and make predictiion about the 
testing data accroding to the results of the cross validation

data_dict contains
    data_dict['percentages']
    data_dict['all_training_distance_matrix']
    data_dict['all_training_design_matrix']
    data_dict['all_training_observed_values']
    data_dict['list_centers']
    data_dict['list_n_centers']
    data_dict['training_set']
    data_dict['testing_set']
    data_dict['positive_training_set']
    data_dict['negative_training_set']
    data_dict['positive_testing_set']
    data_dict['negative_testing_set']
    data_dict['cross_number']


### REMARK: the row number of coefficients and the column number of the design_matrix are the same
    
******************Explanation about the above values***********************
percentages:
    a list of percentages, which are the hyperparameters to be selected by the cross validation
    
all_training_design_matrix:
    is the design matrix with the training set as the rows and the same training set as the centers

all_training_observed_values:
    the obseved values of the training_set, which is the concatenation of the positive training set 
    and the negative training training set.
list_centers:
    a list of selected centers corresponding to the list of percentages. The element of this list
    are integers, which give the position indices in the training set.

list_coeff:
    a list of coefficients of the centers, corresponding to the list_centers

training_set:
    the training set, it contains both the positive training set and the negative training set. The positive
    training set is ahead of the negative training set in the list.

testing_set:
    the testing set, it contains both the positive testing set and the negative testing set and the positive 
    testing set is ahead of the negative testing set in the list.
    
list_AUC:
    a list of areas under different recall precisions curves corresponding to list_centers
    
centers:
    are the centers in a state of constant changing as the selection of the centers progresses.
    It begins with the non_redundant_training set
coeff:
    a list of coefficient, it changes as the process of center selection is going on.
reg:
    regularization coefficient, it changes as the process of the center selection is going on.
    
'''
#########################################################################
def Cross_validation(data_dict, train_para):
    
    # Take out the values
    all_training_distance_matrix = data_dict['all_training_distance_matrix']
    all_training_design_matrix = data_dict['all_training_design_matrix']
    all_training_observed_values = data_dict['all_training_observed_values']# pay attention to this one, it is used in the one step reduce
    positive_training_set = data_dict['positive_training_set']
    negative_training_set = data_dict['negative_training_set']
    training_set = data_dict['training_set']
    cross_number = data_dict['cross_number']
    percentages = data_dict['percentages']
    
    # Initiate list_area_under_rp_curve 
    list_AUC = [0*i for i in range(len(percentages))]
        
    # Generate the cross indices
    cross_train_indices, cross_test_indices = \
    Raised_cross_indices(positive_training_set, negative_training_set, cross_number)    
    
    # Go into the cross
    for i_cross in range(cross_number):
        print('cross   ' + str(i_cross))
        # Get the test indices and the training indices
        test_indices = cross_test_indices[i_cross]
        train_indices= cross_train_indices[i_cross]
        
        # Get the beginning center indices
        centers_init = []; non_redundant = []
        for indx in train_indices:
            if [training_set[indx][0],training_set[indx][1]] not in non_redundant:
                non_redundant.append([training_set[indx][0],training_set[indx][1]])
                centers_init.append(indx)            
        
        # Load the list_n_centers
        train_para['list_n_centers'] = []
        for per in percentages:
            n_center = math.floor(len(centers_init)*per)
            if n_center >= 1:
                train_para['list_n_centers'].append(n_center)
        
        # Load the list_centers
        distance_matrix_center = all_training_distance_matrix[centers_init,:][:, centers_init]
        train_para['list_centers'] =\
        Coverage_reduce_centers(distance_matrix_center,centers_init,train_para['list_n_centers'])
        
        # Get the observed values
        observed_values_train = all_training_observed_values[train_indices]
        observed_values_train = np.reshape(observed_values_train, (-1,1))
#        parameter['observed_values'] = observed_values 
        observed_values_test = all_training_observed_values[test_indices]
        observed_values_test = np.reshape(observed_values_test, (-1,1))

        # Train the model at different centers
        list_centers = train_para['list_centers']
        for i in range(len(list_centers)):
            centers = list_centers[i]
            print('length of centers d%',len(centers))                       
            
#           # Get the testing distance matrix and the training distance matrix
            design_matrix_test = all_training_design_matrix[test_indices,:][:,centers]
            design_matrix_train = all_training_design_matrix[train_indices,:][:,centers]
            # Attach the constant column
            design_matrix_test = np.hstack((design_matrix_test, np.ones((len(test_indices), 1))))
            design_matrix_train = np.hstack((design_matrix_train, np.ones((len(train_indices), 1))))
#            centers.remove(-1)
            
            # initiate the coeff
#            coeff = np.ones((2*len(centers)+1, 1))
            coeff = np.zeros((len(centers)+1, 1))
            train_para['coefficients'] = coeff
            train_para['design_matrix'] = design_matrix_train
            train_para['observed'] = observed_values_train

#            train_para['reg']
#            train_para['loss_type']  Those values should be given outside
#            train_para['method']
            
            # Pay attention to the termination condition
            '''*****************************************'''
#            termination_grad_norm = 1E-4*len(centers)
            termination = 1E-4*len(centers)
            train_para, loss = Train(train_para, rho=0.8, c = 1e-4, termination = termination)
            
            coeff = np.reshape(train_para['coefficients'], (-1,1)) 
            print('length of coefficinets ', np.shape(coeff))
            linear_coeff = train_para['coefficients'][:len(centers)+1, 0]
#            radius_coeff = parameter['coeff'][len(centers)+1:, 0]

            # Make prediction
#            radius_coeff = np.ones((len(centers), 1))
#            radius_coeff = np.reshape(radius_coeff, (-1,1))
            pred = design_matrix_test.dot(linear_coeff)            
            # Calculate the AUC
            AUC, TPR, FPR = Calculate_AUC(pred, observed_values_test)
            print('AUC:  ',AUC)
            AUC /= cross_number            

            list_AUC[i] += AUC
            
    train_para['list_AUC'] = list_AUC  
#################################################################################
'''
pred: a list of predicted values
observed_values: a list of observed values
Those elements in those two list should correspond with each other
'''
def Calculate_AUC(pred, observed_values):
    # Calculate the number of positive observations
    p_s = 0 ; n_s = 0
    for ob in observed_values:
        if ob == 0 :
            n_s += 1
        elif ob > 0:
            p_s += 1
            
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
        elif match[1] == 0:
            n_negative += 1
            
        TPR.append(n_positive/p_s)
        FPR.append(n_negative/n_s)
    
    # Calculate the AUC
    AUC = 0
    for i in range(1, len(TPR)):
        AUC += TPR[i] * (FPR[i] - FPR[i-1])
    
    return AUC, TPR, FPR

################################################################################
def Observed_values(training_set, binary):
    observed_values = []
    if binary:
        for train in training_set:
            if train[2] >0:
                observed_values.append(1)
            else:
                observed_values.append(train[2])
        
        observed_values = np.array(observed_values)
    
    if not binary:
        # we compress the values by the range of the observed values
        sup = 0
        inf = 10000
        for train in training_set:
            if train[2] > sup:
                sup = train[2]
            if train[2] < inf and train[2]>0:
                inf = train[2]
        rang = sup - inf
        # scale the observed values
        for train in training_set:
            if train[2] >0:
                observed_values.append(train[2]/rang)
            else:
                observed_values.append(train[2])
        
        observed_values = np.array(observed_values)
        
    return observed_values
    


######################################################################################
def Distance_design(positive_d, negative_d, basis_function):
    distance_design ={}
    for i in range(1, 5):
        for j in range(1,5):
            p_name = 'training_'+str(i)+'_'+str(j)+'_0_0_1_2_1perchain'
            os.chdir(positive_d)
            with open(p_name, 'r') as f:
                positive_training_set = json.load(f)
            n_name = p_name+'_negative'
            os.chdir(negative_d)
            with open(n_name, 'r') as f:
                negative_training_set = json.load(f)
            
            '''*******The following block is to calculate the observed values 
            ***********according whether the observed should be binary or not'''
            # Calculate 
            training_set = copy.deepcopy(positive_training_set)
            training_set.extend(negative_training_set)   
            # Get the observed values
#            all_training_observed_values = Observed_values(training_set, binary)
            
            print('Working on: ' +str(i)+'_'+str(j)+'_0_0_1_2_1perchain')
            all_training_distance_matrix = Distance_matrix(training_set, training_set, square = True)
            
            radius_coeff = np.ones((all_training_distance_matrix.shape[1], 1))
            all_training_design_matrix = Design_matrix(all_training_distance_matrix,\
                                                       radius_coeff, basis_function)
            
            # load the results
            distance_design[str(i)+'_'+str(j)] = {}
            distance_design[str(i)+'_'+str(j)]['positive_training_set']=positive_training_set
            distance_design[str(i)+'_'+str(j)]['negative_training_set']=negative_training_set
            distance_design[str(i)+'_'+str(j)]['training_set']=training_set
            distance_design[str(i)+'_'+str(j)]['all_training_distance_matrix']=all_training_distance_matrix
            distance_design[str(i)+'_'+str(j)]['all_training_design_matrix']=all_training_design_matrix
            
    return distance_design
###################################################################################
def Batch_cross_validation(distance_design, binary):
      
#                
    cross_number = 5      
    percentages = [0.8, 0.4, 0.2, 0.1, 0.05, 0.02]
#    reg_list = [0, 0.01, 0.05, 0.1, 0.5]
    reg_list = [0.0001, 0.001, 0.01, 0.1]
    
    # The final result will be stored in cross_coverage_RBFN
    cross_coverage_RBFN = {}
    cross_coverage_RBFN['percentages'] = percentages
    cross_coverage_RBFN['cross_number'] = cross_number
    cross_coverage_RBFN['reg_list'] = reg_list    

    for i in range(1, 5):
        for j in range(1,5):
            k = str(i)+'_'+str(j)

            # Load up the parameter
            data_dict = {}
            data_dict['all_training_distance_matrix'] = distance_design[k]['all_training_distance_matrix']
            data_dict['all_training_design_matrix'] = distance_design[k]['all_training_design_matrix']
#            data_dict['all_training_observed_values'] = all_training_observed_values 
            data_dict['positive_training_set'] = distance_design[k]['positive_training_set'] 
            data_dict['negative_training_set'] = distance_design[k]['negative_training_set'] 
            data_dict['training_set'] = distance_design[k]['training_set'] 
            
            data_dict['all_training_observed_values'] =\
                Observed_values(data_dict['training_set'], binary)
            
            data_dict['cross_number'] = cross_number 
            data_dict['percentages'] = percentages 
            
            for reg in reg_list:
                
                train_para ={}
                train_para['reg'] = reg
                train_para['method'] = 'BFGS'
                if binary:
                    train_para['loss_type'] = 'Sigmoid'
                else:
                    train_para['loss_type'] = 'SumSquares'
                
                # Cross_validate
                Cross_validation(data_dict, train_para)               
                # Load the results to cross_coverage_RBFN
                key = str(i)+'_'+str(j)+'_0_0_1_2_1perchain'
                AUC_array = np.array(train_para['list_AUC'])
                AUC_array = np.reshape(AUC_array, (1, -1))
                if key not in cross_coverage_RBFN:
                    cross_coverage_RBFN[key] = AUC_array
                else:                    
                    cross_coverage_RBFN[key] = np.vstack((cross_coverage_RBFN[key], AUC_array))
        
    return cross_coverage_RBFN
'''##################################################################################'''
def Best_parameter(saving_d):
    # The best parameter will be stored in a dictionary
    best_parameter = {}
    # Select the best gap_score and the extened gap score according to the sum of all tests
    os.chdir(saving_d)
    for bi in ['binary', 'numerical']:
        prefix = 'cross_'+bi+'_Gaussian_'
        best_parameter[prefix] = {}
        for gap_score in range(-5, 0):
            for extended_gap_score in range(gap_score, 0):
                file = prefix + str(gap_score)+'_'+str(extended_gap_score)
                with open(file, 'r') as f:
                    results = json.load(f)
                #Calculate the sum of the AUC
                AUC_sum = 0; AUC_largest = 0 
                for i in range(1, 5):
                    for j in range(1, 5):
                        key = str(i)+'_'+str(j)+'_0_0_1_2_1perchain'
                        AUC_lists = results[key]
                        # Take the sum 
                        AUC_sum += np.sum(np.array(AUC_lists))
                if AUC_sum > AUC_largest:
                    AUC_largest = AUC_sum
                    best_gap_ext = [gap_score, extended_gap_score]
        # Load the results
        best_parameter[prefix] ['best_gap_ext'] = best_gap_ext
              
        # Select the best reg and the precentages
        use_file = prefix + str(best_gap_ext[0])+'_'+str(best_gap_ext[1])
        with open(use_file, 'r') as f:
            use_results = json.load(f)
        for i in range(1, 5):
            for j in range(1, 5):
                key = str(i)+'_'+str(j)+'_0_0_1_2_1perchain'
                reg_list = use_results['reg_list']
                percentages = use_results['percentages']
                use_AUC_lists = use_results[key]
                use_max_auc = np.max(np.array(use_AUC_lists))
                best_ind =np.where(np.array(use_AUC_lists) == use_max_auc)
                best_ind = [best_ind[0][0], best_ind[1][0]]
                best_reg_per_auc = [reg_list[best_ind[0]], percentages[best_ind[1]], use_max_auc]
                
                # Save the parameter
                best_parameter[prefix][key]=best_reg_per_auc
            
    return best_parameter

'''##########################################################################'''
'''
Run the main and Save the results
'''
class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)
###################################################################################
#def main(negative_d, positive_d, saving_d):
#    for basis_function in ['Gaussian']:
#        for gap_score in range(-5, 0):
#            for extended_gap_score in range(gap_score, 0):
#            
#                aligner.open_gap_score = gap_score
#                aligner.extend_gap_score = extended_gap_score
#                
#                distance_design = Distance_design(positive_d, negative_d, basis_function)
#                
#                for binary in [True, False]:
#                
#                    cross_coverage_RBFN_binary = \
#                    Batch_cross_validation(distance_design, binary)
#                    
#                    os.chdir(saving_d)
#                    suffix = '_'+str(gap_score)+'_'+str(extended_gap_score)
#                    if binary:
#                        name = 'cross_binary_'+basis_function+suffix
#                    if not binary:
#                        name = 'cross_numerical_'+basis_function+suffix
#                    with open(name, 'w') as f:
#                        json.dump(cross_coverage_RBFN_binary, f, cls=NumpyEncoder)
#                        
#    # Find the best hyperparameter and save
#    best_parameter = Best_parameter(saving_d)
#    os.chdir(saving_d)
#    with open('best_parameter', 'w') as f:
#        json.dump(best_parameter, f)
        
'''***************************************************************************'''
'''***************************************************************************'''
'''***************************************************************************'''
'''***************************************************************************'''

##############################################################################
if __name__ == '__name__':
    positive_d = '/home/leo/Documents/Database/Data_Code_Publish/Cores/Positive_cores'
    negative_d = '/home/leo/Documents/Database/Data_Code_Publish/Cores/Negative_cores/Sample_0'
    saving_d = '/home/leo/Documents/Database/Data_Code_Publish/Codes/Cross_validation'
    
    for basis_function in ['Gaussian']:
        for gap_score in range(-5, 0):
            for extended_gap_score in range(gap_score, 0):
            
                aligner.open_gap_score = gap_score
                aligner.extend_gap_score = extended_gap_score
                
                distance_design = Distance_design(positive_d, negative_d, basis_function)
                
                for binary in [True, False]:
                
                    cross_coverage_RBFN_binary = \
                    Batch_cross_validation(distance_design, binary)
                    
                    os.chdir(saving_d)
                    suffix = '_'+str(gap_score)+'_'+str(extended_gap_score)
                    if binary:
                        name = 'cross_binary_'+basis_function+suffix
                    if not binary:
                        name = 'cross_numerical_'+basis_function+suffix
                    with open(name, 'w') as f:
                        json.dump(cross_coverage_RBFN_binary, f, cls=NumpyEncoder)
                        
    # Find the best hyperparameter and save
    best_parameter = Best_parameter(saving_d)
    os.chdir(saving_d)
    with open('best_parameter', 'w') as f:
        json.dump(best_parameter, f)
##############################################################################
#
#def Sub_test_AUC(testing_set, centers_selected, coeff, binary, basis_function):
#    
#    linear_coeff = coeff[:len(centers_selected)+1, 0]
#    radius_coeff = coeff[len(centers_selected)+1:, 0]
#    
#    linear_coeff = np.reshape(linear_coeff, (-1, 1))
#    radius_coeff = np.reshape(radius_coeff, (-1, 1))
#    
#    testing_distance_matrix = Distance_matrix(testing_set, centers_selected, square = False)
#    testing_design_matrix = Design_matrix(testing_distance_matrix, radius_coeff, basis_function)
#    row_n = np.shape(testing_design_matrix)[0]
#    testing_design_matrix = np.hstack((testing_design_matrix, np.ones((row_n, 1))))
#    
#    pred = testing_design_matrix.dot(linear_coeff)     
#
#    testing_observed_values = Observed_values(testing_set, binary)     
#    # Calculate the AUC
#    AUC, TPR, FPR = Calculate_AUC(pred, testing_observed_values)
#    print('AUC:  ',AUC)         
#    
#    return AUC, TPR, FPR
##########################################################################
#'''
#From the cross_coverage_RBFN, we see that more the number of the centers, the better
#the prediction. Thus we use all the non_redundent_training set as the centers
#
#Test_coverage_RBFN:
#    This function is to calculate the precision recall rate for the positive_testing_set
#    and the negative_testing_set. There are ten different sets of negative_testing_set. We
#    can calculate the separately
#Input:
#    best_percentage:
#        The percentage which performs the best in the cross_coverage_RBFN
#Output:
#    test_coverage_RBFN_results:
#        a list of dictionaries, each dictionary contains the following values
#        coeff:
#            a list of coefficients
#        non_redundant_training_set:
#            a sub set of training_set with the redundant removed
#        area_list:
#            a list of areas, for each negative_testing_set
#        recall_list:
#            a list of recalls, each recall is a list, corresponding to one 
#            set of negative_testing_set.
#        precision_list:
#            a list of precisions, each precision is a list, corresponding to one
#            set of negative_testing_set.
#'''
#def Test_AUC( best_hyperparameter,binary = True,basis_function='Gaussian'):
#    test_results={}
#    # Read the training_data and the testing data
#   
#    negative_d = '/home/leo/Documents/Database/Pipeline_New/Complexes/Negative_cores/Sample_0'
#    positive_d = '/home/leo/Documents/Database/Pipeline_New/Complexes/Cores'
#    n_directory = '/home/leo/Documents/Database/Pipeline_New/Complexes/Negative_cores/Sample_'
#      
#    for i in range(1, 4):
#        for j in range(1, 4):
#            # set the empty container
#            key = str(i)+'_'+str(j)+'_0_0_1_2_1perchain'
#            test_results[key] = {}
#            test_results[key]['AUC_average'] = 0 
#            test_results[key]['TPR_list'] = []
#            test_results[key]['FPR_list'] = []
#            
#            p_train = 'training_'+str(i)+'_'+str(j)+'_0_0_1_2_1perchain'
#            p_test = 'testing_'+str(i)+'_'+str(j)+'_0_0_1_2_1perchain'
#            os.chdir(positive_d)
#            with open(p_train, 'r') as f:
#                positive_training_set = json.load(f)
#                
#            n_traing = p_train+'_negative'
#            os.chdir(negative_d)
#            with open(n_traing, 'r') as f:
#                negative_training_set = json.load(f)
#                
#            training_set = copy.deepcopy(positive_training_set)
#            training_set.extend(negative_training_set)
#            
#            # Get the best hyperparameters
#            best = best_hyperparameter[str(i)+'_'+str(j)]
#            percentage = best[1]
#            binary = False
#            basis_function = 'Gaussian'
#            reg = 0
#            
#            
#            
#            print('Working on '+ key)  
#            parameter, centers_selected =\
#                Train_use_best_hyperparameter(training_set,percentage, reg,binary, basis_function)
#                
#            coeff = np.reshape(parameter['coeff'], (-1,1)) 
#            test_results[key]['coeff'] = coeff
#            test_results[key]['centers_selected']=centers_selected
#            print('length of coefficinets ', np.shape(coeff))
#            
#            # Load the positive testing set
#            os.chdir(positive_d)
#            with open(p_test, 'r') as f:
#                positive_testing_set = json.load(f)
#            for n in range(10):
#                testing_set = copy.deepcopy(positive_testing_set)
#                os.chdir(n_directory+str(n))
#                with open(p_test+'_negative', 'r') as f:
#                    negative_testing_set = json.load(f)
#                testing_set.extend(negative_testing_set)
#                AUC, TPR, FPR = Sub_test_AUC(testing_set, \
#                                             centers_selected, coeff, binary, basis_function)
#                test_results[key]['AUC_average'] += AUC/10
#                test_results[key]['TPR_list'].append(TPR)
#                test_results[key]['FPR_list'].append(FPR)
#            
#    return test_results
#
#
##test_results.keys()
##test_results['1_1_0_0_1_2_1perchain']['AUC_average']
##test_results['1_1_0_0_1_2_1perchain']['TPR_list']
##len(test_results['1_1_0_0_1_2_1perchain']['coeff'])
##len(test_results['1_1_0_0_1_2_1perchain']['centers_selected'])
################################################################################

###############################################################################
#'''
#Reverse_test:
#    A function to test the testing set with the Ab and Ag sequences switched
#Input:
#    test_results
#Output:
#    reverse_test_coverage_RBFN_results
#    ['recalls_list']
#    [precidions_list]
#    [area_list]
#'''
#def Discrimation_test(test_results):
#    reverse_test_results={}
#    # Read the training_data and the testing data
#    positive_d = '/home/leo/Documents/Database/Pipeline_New/Complexes/Cores'
#    for i in range(1, 4):
#        for j in range(1, 4):
#            # set the empty container
#            key = str(i)+'_'+str(j)+'_0_0_1_2_1perchain' 
#            reverse_test_results[key] = {}
#
#            p_test = 'testing_'+str(i)+'_'+str(j)+'_0_0_1_2_1perchain'
#            os.chdir(positive_d)
#            with open(p_test, 'r') as f:
#                positive_testing_set = json.load(f)  
#            
#            os.chdir(positive_d)
#            n_test = 'testing_'+str(j)+'_'+str(i)+'_0_0_1_2_1perchain'
#            with open(n_test, 'r') as f:
#                negative_testing = json.load(f)
#
#            negative_testing_set = []
#            for parepi_negative in negative_testing:
#                negative_testing_set.append([parepi_negative[1], parepi_negative[0],\
#                                             -1, parepi_negative[3]])                    
#            
#            testing_set = copy.deepcopy(positive_testing_set)
#            testing_set.extend(negative_testing_set)
#            
#            coeff = test_results[key]['coeff']
#            coeff = np.reshape(np.array(coeff), (-1,1))
#            centers_selected = test_results[key]['centers_selected']
#            AUC, TPR, FPR = Sub_test_AUC(testing_set,\
#                        centers_selected, coeff, binary=False, basis_function='Gaussian')
#                
#            reverse_test_results[key]['AUC'] = AUC
#            reverse_test_results[key]['TPR'] = TPR
#            reverse_test_results[key]['FPR'] = FPR
#            
#    return reverse_test_results
########################################################################
#def All_tests():
#    best = Best_parameter()
#    test_results = Test_AUC(best,binary = False,basis_function='Gaussian')
#    reverse_test_results = Discrimation_test(test_results)
#    os.chdir('/home/leo/Documents/Database/Pipeline_New/Complexes/Results')
#    with open('test_results', 'w') as f:
#        json.dump(test_results, f, cls=NumpyEncoder)
#    with open('reverse_test_results', 'w') as f:
#        json.dump(reverse_test_results, f, cls=NumpyEncoder)
#    
#    
##All_tests()        
#
#        
##    reverse_test_results = Discrimation_test()
##    os.chdir('/home/leo/Documents/Database/Pipeline_New/Results')
##    with open('reverse_test_results', 'w') as f:
##        json.dump(reverse_test_results, f)
########################################################################
##All_tests()




