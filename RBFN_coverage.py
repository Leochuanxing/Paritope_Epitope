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
os.chdir('/home/leo/Documents/Database/Pipeline_New/Codes')
from AlignmentConstraint import To_seq

#############################################################
# Define  distances
from Bio import Align
from Bio.SubsMat.MatrixInfo import blosum62
aligner = Align.PairwiseAligner()
aligner.substitution_matrix = blosum62
aligner.open_gap_score = -5
aligner.extend_gap_score = -1
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
    return math.exp(-distance**2/radius**2)
    
def Mrakov(distance, radius):
    return math.exp(-distance/radius)
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
'''
def Design_matrix(distance_matrix, basis_function = 'Markov', radius = 1):
    nrow = np.shape(distance_matrix)[0]
    ncol = np.shape(distance_matrix)[1]
    
    design_matrix = np.zeros_like(distance_matrix)
    for i in range(nrow):
        for j in range(ncol):
            if basis_function == 'Gaussian':
                design_matrix[i, j] = Gaussian(distance_matrix[i, j], radius)
            elif basis_function == 'Markov':
                design_matrix[i, j] = Mrakov(distance_matrix[i, j], radius)
            elif basis_function == 'Inverse_Multi_Quadric':
                design_matrix[i, j] = Inverse_Multi_Quadric(distance_matrix[i, j], c=1, beta=2)  
                
    # Add a constant term to the design matrix
    design_matrix = np.hstack((design_matrix, np.ones((nrow,1))))
           
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

def Coverage_reduce_centers(training_training_distance_matrix, training_set,  n_centers_list):
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
Use a small distance matrix to check the above function, which is very important.
'''
#os.chdir('/home/leo/Documents/Database/Pipeline_New/Cores')
#with open('training_2_2_1_1_1_2_1perchain', 'r') as f:
#    training_set = json.load(f)
#len(training_set)
#small_data = training_set[6:12]
#small_distance_matrix = Distance_matrix(small_data, small_data, square = True)
#small_distance_matrix = np.round(small_distance_matrix, 3)
#small_distance_matrix
#n_centers_list = [1,2,3,4,5,6]
#small_training_set = ['a', 'b', 'c', 'd', 'e', 'f']
#centers_list= Coverage_reduce_centers(small_distance_matrix,small_training_set, n_centers_list)
#centers_list

            
##############################################################################
'''
Remove the duplicates before removing the centers
'''
def Remove_duplicates(training_set):
    centers = []
    pure_parepi = []
    non_redundent_training_set = []
    for i in range(len(training_set)):
        parepi = training_set[i]
        if [parepi[0], parepi[1]] not in pure_parepi:            
            pure_parepi.append([parepi[0], parepi[1]])
            centers.append(i)
            non_redundent_training_set.append(training_set[i])
    return centers, non_redundent_training_set
##############################################################################
'''
Loss:
    To return the loss, the gradient, and the parameters
Input:
        The return value of Design_matrix, in the shape of (m,n)
    observed_values:
        A vector gives the observed valuse, in the shape of (m,1)
    parameter:
        A dictionary, 
        parameter['coeff'] contains the vector of the coefficients, in the shape of (n,1)
        parameter['reg'] contains the vector of the regularization coefficients, in the shape of (n,1)
        Here we had better design this loss function so that it can be used in the 
        ridge selection of the centers as well.        
Output:
    loss:
        A real number
    gradient:
        a dictionary
        gradient['coeff'] is the gradient of the parameter['coeff']
        gradient['reg']  is the gradient of the parameter['reg']
'''
def Loss(design_matrix, observed_values, parameter):
    # Unpack the dictionary 
    coeff = parameter['coeff']
    reg = parameter['reg']
    
    coeff_square = coeff*coeff
    diff = design_matrix.dot(coeff) - observed_values
    loss = (diff.T).dot(diff)
    loss += (coeff_square.T).dot(reg)
    grad_coeff = 2*(design_matrix.T).dot(diff) + 2 * coeff * reg
    # Pack up the results
    gradient = {}
    gradient['coeff'] = grad_coeff
    
    return loss, gradient
####################################################################################

'''
Train_RBFN_BFGS:
    A function to train the RBFN by using the BFGS method
Input:
    design_matrix:
        The return value of Design_matrix, in the shape of (m,n)
    observed_values:
        A vector gives the observed valuse, in the shape of (m,1)
    rho:
        The contraction factor in Ramijo backtracking
    c:
        the Goldstein coefficient c
    termination:
        The termination condition with the norm of the gradient < termination
Output:
    parameter:
        A dictionary contains
        parameter['coeff'], the coefficient after training
        parameter['reg'], the regularization coefficient    
'''


def Train_RBFN_BFGS(design_matrix, observed_values, rho=0.9, c = 1e-3, termination = 1e-2,\
                    parameter_inheritance = False, parameter=None):
    
    nrow = np.shape(design_matrix)[0]
    ncol = np.shape(design_matrix)[1]
        
    # Give the initial Hessian H. The variables are the coeff and reg
    H = np.eye(ncol)/(10*nrow)
    # Check if it inherit the parameter from somewhere else.
    if not parameter_inheritance :
        # Set the starting point
        coeff = np.zeros((ncol,1))
        parameter = {}
        parameter['coeff'] = coeff
        #The reg should not be negative. It is better that reg > delta, a small positive number
        reg = np.ones((ncol,1)) 
        parameter['reg'] = reg

    # BFGS algorithm
    loss, gradient = Loss(design_matrix, observed_values, parameter)
    grad_coeff = gradient['coeff']
    grad_coeff = np.reshape(grad_coeff, (-1, 1))
    ternination_square = termination**2
    grad_square = ternination_square + 1
#    grad_square = (grad_coeff.T).dot(grad_coeff)
    while grad_square >= ternination_square:        
        p = - H.dot(grad_coeff)        
        # Find the next coeff
        parameter_new = {}
        parameter_new['coeff'] = p + parameter['coeff']
        parameter_new['reg'] = parameter['reg']
        
        new_loss, new_gradient = Loss(design_matrix, observed_values, parameter_new)
        # Ramijo Back-tracking
        while new_loss > loss + c * (grad_coeff.T).dot(p):
            p *= rho
            parameter_new['coeff'] = p + parameter['coeff']            
            new_loss, new_gradient = Loss(design_matrix, observed_values, parameter_new)
        
        # update H
        s = p
        new_grad = new_gradient['coeff']
        y = new_grad - grad_coeff
        r = (y.T).dot(s)
        if r != 0:
            r = 1/r
            I = np.eye(ncol)
            H = (I - r*s.dot(y.T)).dot(H).dot(I - r*y.dot(s.T)) + r*s.dot(s.T)# Can be accelerate
        else:
            H = I
        # Update loss, grad, grad_square and paramter
        loss = new_loss
        grad_coeff = new_grad
        parameter['coeff'] = parameter_new['coeff']
        grad_square = (grad_coeff.T).dot(grad_coeff)
        print('loss  ', loss, '    ','grad_square   ', grad_square)
    return parameter, loss
###############################################################################
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
    cross_training_set:
        a list contains 'cross_number' of lists of positions of the training_set
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

parameter contains
    parameter['percentages']
    parameter['all_training_distance_matrix']
    parameter['all_training_design_matrix']
    parameter['all_training_observed_values']
    parameter['best_coeff']
    parameter['best_reg']
    parameter['best_centers']
    parameter['list_centers']
    parameter['list_n_centers']
    parameter['list_coeff']
    parameter['list_reg']
    parameter['training_set']
    parameter['testing_set']
    parameter['positive_training_set']
    parameter['negative_training_set']
    parameter['positive_testing_set']
    parameter['negative_testing_set']
    parameter['cross_number']
    parameter['list_area_under_rp_curve']
    parameter['best_recall_precision_list']
    parameter['design_matrix']
    parameter['observed_values']
    parameter['centers']
    parameter['coeff']
    parameter['reg']
    
******************Explanation about the above values***********************
percentages:
    a list of percentages, which are the hyperparameters to be selected by the cross validation
    
all_training_design_matrix:
    is the design matrix with the training set as the rows and the same training set as the centers

all_training_observed_values:
    the obseved values of the training_set, which is the concatenation of the positive training set 
    and the negative training training set.
    
observed_values:
    a list of observed values used to train the model and do center selection. They are corresponding 
    to the rows of the design matrix
    
best_coeff:
    After the best parameter is selected, the model is trained again to ge this best_coeff
    
best_reg:
    the best regularization coefficienct, in this model, the reg is constant 1s.
    
best_centers:
    gives the best centers selected. Those centers are concrete instead of the position indices.

list_centers:
    a list of selected centers corresponding to the list of percentages. The element of this list
    are integers, which give the position indices in the training set.

list_coeff:
    a list of coefficients of the centers, corresponding to the list_centers

list_reg:
    a list of regulariztion coefficients, corresponding to the list_centers
    
training_set:
    the training set, it contains both the positive training set and the negative training set. The positive
    training set is ahead of the negative training set in the list.

testing_set:
    the testing set, it contains both the positive testing set and the negative testing set and the positive 
    testing set is ahead of the negative testing set in the list.
    
list_area_under_rp_curve:
    a list of areas under different recall precisions curves corresponding to list_centers

best_recall_precision_list:
    a list contains two lists, one list are the recall rates, the other list are the precisions
    those two lists correspond to each other, and they are calculated from the best_coeff
design_matrix:
    this design matrix is changing as the selection of the centers progresses. 
centers:
    are the centers in a state of constant changing as the selection of the centers progresses.
    It begins with the non_redundant_training set
coeff:
    a list of coefficient, it changes as the process of center selection is going on.
reg:
    regularization coefficient, it changes as the process of the center selection is going on.
    
'''
#########################################################################
def Cross_validation(parameter):

    # Take out the values
    all_training_distance_matrix = parameter['all_training_distance_matrix ']
    all_training_design_matrix = parameter['all_training_design_matrix']
    all_training_observed_values = parameter['all_training_observed_values']# pay attention to this one, it is used in the one step reduce
    positive_training_set = parameter['positive_training_set']
    negative_training_set = parameter['negative_training_set']
    training_set = parameter['training_set']
    cross_number = parameter['cross_number']
    percentages = parameter['percentages']
    
    # Initiate list_area_under_rp_curve 
    list_area_under_rp_curve = [0*i for i in range(len(percentages))]
        
    # Generate the cross indices
    cross_train_indices, cross_test_indices = \
    Raised_cross_indices(positive_training_set, negative_training_set, cross_number)
    
    
    # Go into the cross
    for i_cross in range(cross_number):
        print('cross   ' + str(i_cross))
        # Get the row_indices and the col_indices
        row_indices = cross_test_indices[i_cross]
        col_indices= cross_train_indices[i_cross]
        
        # Get the cross training set
        cross_training_set = []
        for indx in col_indices:
            cross_training_set.append(training_set[indx])
                
        # Load the beginning centers
        centers, non_redundent_training_set = Remove_duplicates(cross_training_set)
        parameter['centers']=centers
        
        # Load the list_n_centers
        parameter['list_n_centers'] = []
        for per in percentages:
            n_center = math.floor(len(centers)*per)
            if n_center >= 1:
                parameter['list_n_centers'].append(n_center)
        
        # Load the list_centers
        distance_matrix = all_training_distance_matrix[centers,:][:, centers]
        parameter['list_centers'] =\
        Coverage_reduce_centers(distance_matrix,centers,parameter['list_n_centers'])
        
        # Load the observed values
        observed_values = all_training_observed_values[row_indices]
        observed_values = np.reshape(observed_values, (-1,1))
        parameter['observed_values'] = observed_values        

        # Train the model at different centers
        list_centers = parameter['list_centers']
        for i in range(len(list_centers)):
            centers = list_centers[i]
            print(len(centers))
            
            centers.append(-1)
            design_matrix = all_training_design_matrix[row_indices,:][:,centers]
            centers.remove(-1)
            
            # initiate the coeff
            coeff = np.zeros((len(centers)+1, 1))
            parameter['coeff'] = coeff
            # Initiate the reg
            reg = np.ones((len(coeff),1))
            parameter['reg'] = reg
            
            # Pay attention to the termination condition
            parameter, loss = Train_RBFN_BFGS(design_matrix, observed_values, rho=0.9, c = 1e-3, termination = 1e-2,\
                                             parameter_inheritance = True, parameter=parameter)
            
            coeff = parameter['coeff']
            # Make prediction
            pred = design_matrix.dot(coeff)
            # Match up the observed values with the pred
            match_up = []
            for j in range(len(observed_values)):
                match_up.append([pred[j], observed_values[j]])
            match_up.sort(key = lambda x:x[0], reverse = True)
            
            area = 0
            n_positive = 0
            n_negative = 0
            for match in match_up:
                if match[1] == 1:
                    n_positive += 1
                elif match[1] == -1:
                    n_negative += 1
                area += n_positive/(n_positive+n_negative)
            
            # Take the average
            area/=(len(row_indices) * cross_number)
            # load to the list area under rp curve
            list_area_under_rp_curve[i] += area
            
    parameter['list_area_under_rp_curve'] = list_area_under_rp_curve  
#################################################################################
#'''
#Try the cross validation on one set of data
#'''
#os.chdir('/home/leo/Documents/Database/Pipeline_New/Cores')  
#with open('training_2_2_1_1_1_2_1perchain', 'r') as f:
#    positive_training_set = json.load(f)
#os.chdir('/home/leo/Documents/Database/Pipeline_New/Negative_Cores/Sample_0')
#with open('training_2_2_1_1_1_2_1perchain_negative', 'r') as f:
#    negative_training_set = json.load(f)
#
## Calculate 
#training_set = copy.deepcopy(positive_training_set)
#training_set.extend(negative_training_set)
#
#all_training_observed_values = []
#for train in training_set:
#    if train[2] >0:
#        all_training_observed_values.append(1)
#    else:
#        all_training_observed_values.append(train[2])
#all_training_observed_values = np.array(all_training_observed_values)
#
#print('Calculating the distance matrix')
#all_training_distance_matrix = Distance_matrix(training_set, training_set, square = True)
#
#all_training_design_matrix = Design_matrix(all_training_distance_matrix)
#
#cross_number = 6
#
#percentages = [1, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001]
## Load up the parameter
#parameter = {}
#parameter['all_training_distance_matrix '] = all_training_distance_matrix 
#parameter['all_training_design_matrix'] = all_training_design_matrix
#parameter['all_training_observed_values'] = all_training_observed_values 
#parameter['positive_training_set'] = positive_training_set 
#parameter['negative_training_set'] = negative_training_set 
#parameter['training_set'] = training_set 
#parameter['cross_number'] = cross_number 
#parameter['percentages'] = percentages 
#
## Cross_validate
#Cross_validation(parameter)
#
#parameter['list_area_under_rp_curve']
#
#with open('/home/leo/Documents/Database/Pipeline_New/Results/cross_result_RBFN', 'r') as f:
#    cross_RBFN = json.load(f)
#cross_RBFN['2_2_1_1_1_2_1perchain']
###################################################################################
def Batch_cross_validation():
    # The final result will be stored in cross_coverage_RBFN
    cross_coverage_RBFN = {}
    negative_d = '/home/leo/Documents/Database/Pipeline_New/Negative_Cores/Sample_0'
    positive_d = '/home/leo/Documents/Database/Pipeline_New/Cores'
    for i in range(1, 5):
        for j in range(1,5):
            p_name = 'training_'+str(i)+'_'+str(j)+'_1_1_1_2_1perchain'
            os.chdir(positive_d)
            with open(p_name, 'r') as f:
                positive_training_set = json.load(f)
            n_name = p_name+'_negative'
            os.chdir(negative_d)
            with open(n_name, 'r') as f:
                negative_training_set = json.load(f)
                
            # Calculate 
            training_set = copy.deepcopy(positive_training_set)
            training_set.extend(negative_training_set)
            
            all_training_observed_values = []
            for train in training_set:
                if train[2] >0:
                    all_training_observed_values.append(1)
                else:
                    all_training_observed_values.append(train[2])
            all_training_observed_values = np.array(all_training_observed_values)
            
            print('Calculating the distance matrix')
            all_training_distance_matrix = Distance_matrix(training_set, training_set, square = True)
            
            all_training_design_matrix = Design_matrix(all_training_distance_matrix)
            
            cross_number = 6
            
            percentages = [1, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001]
            # Load up the parameter
            parameter = {}
            parameter['all_training_distance_matrix '] = all_training_distance_matrix 
            parameter['all_training_design_matrix'] = all_training_design_matrix
            parameter['all_training_observed_values'] = all_training_observed_values 
            parameter['positive_training_set'] = positive_training_set 
            parameter['negative_training_set'] = negative_training_set 
            parameter['training_set'] = training_set 
            parameter['cross_number'] = cross_number 
            parameter['percentages'] = percentages 
            
            # Cross_validate
            Cross_validation(parameter)
            
            # Load the results to cross_coverage_RBFN
            key = str(i)+'_'+str(j)+'_1_1_1_2_1perchain'
            cross_coverage_RBFN[key] = parameter['list_area_under_rp_curve']
        
    return cross_coverage_RBFN
                
######################################################################################
'''
Run the main and Save the results
'''
class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)
###################################################################################



