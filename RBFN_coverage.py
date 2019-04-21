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
    return math.exp(-(distance**2)*(radius**2))
    
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
        Distance_matrix in the shape of (m,n)
    observed_values:
        A vector gives the observed valuse, in the shape of (m,1)
    parameter:
        A dictionary, 
        parameter['coeff'] contains the vector of the coefficients, in the shape of (2n+1,1)
            where the first n+1 elements are the linear coefficients and the last n elements are
            the radius coefficients.
        parameter['reg']: a float, gives the coeffincient of the regularizaion coefficient
    distance_matrix
       
Output:
    loss:
        A real number
    gradient:
        a dictionary
        gradient['coeff'] is the gradient of the parameter['coeff']
'''
def Loss(distance_matrix, observed_values, parameter):
    # Unpack the dictionary 
    reg = parameter['reg']
    n_row, n_col = distance_matrix.shape
    linear_coeff = parameter['coeff'][:n_col+1, 0]
    linear_coeff = np.reshape(linear_coeff, (-1,1))
    radius_coeff = parameter['coeff'][n_col+1:,0]

    radius_coeff_mat = np.tile(np.reshape(radius_coeff,(1, -1)), (n_row, 1))
    D_square = distance_matrix * distance_matrix
    
    radius_coeff = np.reshape(radius_coeff, (-1,1))
    design_matrix_pure = Design_matrix(distance_matrix,radius_coeff, basis_function = 'Gaussian')
    design_matrix = np.hstack((design_matrix_pure, np.ones((n_row,1))))    

    # Calculate the loss
    linear_coeff = np.reshape(linear_coeff, (-1,1))
#    loss = np.sum(linear_coeff*linear_coeff*reg) + np.sum(radius_coeff * radius_coeff*reg)
    loss = np.sum(linear_coeff*linear_coeff*reg)
    diff = design_matrix.dot(linear_coeff) - observed_values
    loss += (diff.T).dot(diff) /n_row

#    loss += np.sum(coeff_square*reg)
    # Calculate the grad of the linear_coeff
    grad_linear_coeff = 2*(design_matrix.T).dot(diff) 
    grad_linear_coeff /= n_row
    grad_linear_coeff += 2 * linear_coeff * reg
    # Calculate the grad of the radius_coeff

    pure_l_c = np.reshape(linear_coeff[:-1,0], (-1,1))
    grad_radius_coeff = np.sum(2*(diff.dot(pure_l_c.T))*design_matrix_pure*D_square*(-2)*radius_coeff_mat, axis=0)
    grad_radius_coeff = np.reshape(grad_radius_coeff, (-1, 1))
    grad_radius_coeff /= n_row
#    grad_radius_coeff += 2*radius_coeff*reg
    # Pack up the results
    gradient = {}
#    gradient['linear_coeff'] = grad_linear_coeff
#    gradient['radius_coeff'] =grad_radius_coeff
    gradient['coeff'] = np.vstack((grad_linear_coeff, grad_radius_coeff))
#    print('Design matrix pure', design_matrix_pure)
#    print('pure_l_c', pure_l_c)
    
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
#a = np.reshape(np.array([1,2,3]), (-1,1))
#b = np.reshape(np.array([4,5,6]), (-1,1))
#c = np.vstack((a, b))
#c.shape
def Train_RBFN_BFGS(distance_matrix, observed_values, rho=0.8, c = 1e-3, termination = 1e-2,\
                    parameter_inheritance = False, parameter=None):
    
    nrow = np.shape(distance_matrix)[0]
    ncol = np.shape(distance_matrix)[1]
        
    # Give the initial Hessian H. The variables are the coeff and reg
    H = np.eye(2*ncol+1)/(10*nrow)
    # Check if it inherit the parameter from somewhere else.
    if not parameter_inheritance :
        # Set the starting point
        parameter = {}
#        parameter['linear_coeff'] = np.zeros((ncol+1,1))
#        parameter['radius_coeff'] = np.ones((ncol,1))
        parameter['coeff'] = np.ones((2*ncol+1,1))
        #The reg should not be negative. It is better that reg > delta, a small positive number
#        parameter['reg'] = reg

    # BFGS algorithm
    loss, gradient = Loss(distance_matrix, observed_values, parameter)
#    grad_linear_coeff = gradient['linear_coeff']
#    grad_linear_coeff = np.reshape(grad_linear_coeff, (-1, 1))
#    grad_radius_coeff = gradient['radius_coeff']
#    grad_radius_coeff = np.reshape(grad_radius_coeff, (-1,1))
    grad_coeff = gradient['coeff']
    ternination_square = termination**2
    grad_square = ternination_square + 1
#    grad_square = (grad_coeff.T).dot(grad_coeff)
    while grad_square >= ternination_square:        
        p = - H.dot(grad_coeff)        
        # Find the next coeff
        parameter_new = {}
        parameter_new['coeff'] = p + parameter['coeff']
        parameter_new['reg'] = parameter['reg']
        
        new_loss, new_gradient = Loss(distance_matrix, observed_values, parameter_new)
        # Ramijo Back-tracking
        while new_loss > loss + c * (grad_coeff.T).dot(p):
            p *= rho
            parameter_new['coeff'] = p + parameter['coeff']            
            new_loss, new_gradient = Loss(distance_matrix, observed_values, parameter_new)
        
        # update H
        s = p
        new_grad = new_gradient['coeff']
        y = new_grad - grad_coeff
        r = (y.T).dot(s)
        I = np.eye(2*ncol+1)
        if r != 0:
            r = 1/r            
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
def Train_RBFN_BP(distance_matrix, observed_values,step_size = 1E-8,
                    parameter_inheritance = False, parameter=None):
    
    nrow, ncol = np.shape(distance_matrix)
    # Check if it inherit the parameter from somewhere else.
    if not parameter_inheritance :
        # Set the starting coefficient
        parameter = {}
        parameter['coeff'] = np.ones((2*ncol+1,1))

    n = 0
    while n <1E7:
        n+=1
        loss, gradient = Loss(distance_matrix, observed_values, parameter)
        parameter['coeff'] -= step_size*gradient['coeff']
        print(loss)
        
        
    return parameter, loss

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
    
list_AUC:
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
def Cross_validation(parameter, method):

    # Take out the values
    all_training_distance_matrix = parameter['all_training_distance_matrix ']
#    all_training_design_matrix = parameter['all_training_design_matrix']
    all_training_observed_values = parameter['all_training_observed_values']# pay attention to this one, it is used in the one step reduce
    positive_training_set = parameter['positive_training_set']
    negative_training_set = parameter['negative_training_set']
    training_set = parameter['training_set']
    cross_number = parameter['cross_number']
    percentages = parameter['percentages']
    
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
        centers = []; non_redundant = []
        for indx in train_indices:
            if [training_set[indx][0],training_set[indx][1]] not in non_redundant:
                non_redundant.append([training_set[indx][0],training_set[indx][1]])
                centers.append(indx)
            
                
        # Load the beginning centers
#        centers_relative, non_redundent_training_set = Remove_duplicates(cross_training_set)
        # the indices in the above center are corresponding to cross_training_set
        # they should be mapped to the original indices
#        centers = []
#        for c in centers_relative:
#            centers.append(train_indices[c])
            
#        parameter['centers']=centers
        
        # Load the list_n_centers
        parameter['list_n_centers'] = []
        for per in percentages:
            n_center = math.floor(len(centers)*per)
            if n_center >= 1:
                parameter['list_n_centers'].append(n_center)
        
        # Load the list_centers
        distance_matrix_center = all_training_distance_matrix[centers,:][:, centers]
        parameter['list_centers'] =\
        Coverage_reduce_centers(distance_matrix_center,centers,parameter['list_n_centers'])
        
        # Get the observed values
        observed_values_train = all_training_observed_values[train_indices]
        observed_values_train = np.reshape(observed_values_train, (-1,1))
#        parameter['observed_values'] = observed_values 
        observed_values_test = all_training_observed_values[test_indices]
        observed_values_test = np.reshape(observed_values_test, (-1,1))

        # Train the model at different centers
        list_centers = parameter['list_centers']
        for i in range(len(list_centers)):
            centers = list_centers[i]
            print('length of centers d%',len(centers))                       
            
#           # Get the testing distance matrix and the training distance matrix
            distance_matrix_testing = all_training_distance_matrix[test_indices,:][:,centers]
            distance_matrix_training = all_training_distance_matrix[train_indices,:][:,centers]
#            centers.remove(-1)
            
            # initiate the coeff
            coeff = np.ones((2*len(centers)+1, 1))
            parameter['coeff'] = coeff
            
            # Pay attention to the termination condition
            '''*****************************************'''
            termination_grad_norm = 1E-4*len(centers)
            if method == 'BFGS':
                parameter, loss = Train_RBFN_BFGS(distance_matrix_training, observed_values_train, rho=0.5, c = 1e-3,\
                                                  termination = termination_grad_norm,\
                                                 parameter_inheritance = True, parameter=parameter)
#            elif method =='BP':
#                parameter, loss = Train_RBFN_BP(distance_matrix_training, observed_values_train, step_size = 1E5,\
#                                                parameter_inheritance=True, parameter=parameter)
            
            coeff = np.reshape(parameter['coeff'], (-1,1)) 
            print('length of coefficinets ', np.shape(coeff))
            linear_coeff = parameter['coeff'][:len(centers)+1, 0]
            radius_coeff = parameter['coeff'][len(centers)+1:, 0]

            # Make prediction
            radius_coeff = np.reshape(radius_coeff, (-1,1))
            design_matrix_test_pure = Design_matrix(distance_matrix_testing, radius_coeff, basis_function='Gaussian')
            row_n = np.shape(design_matrix_test_pure)[0]
            design_matrix_test = np.hstack((design_matrix_test_pure, np.ones((row_n, 1))))
            pred = design_matrix_test.dot(linear_coeff)            
            # Calculate the AUC
            AUC, TPR, FPR = Calculate_AUC(pred, observed_values_test)
            print('AUC:  ',AUC)
            AUC /= cross_number            

            list_AUC[i] += AUC
            
    parameter['list_AUC'] = list_AUC  
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
        if ob < 0 :
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
        elif match[1] == -1:
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
    
###################################################################################
def Batch_cross_validation(binary, method):
    # The final result will be stored in cross_coverage_RBFN
    cross_coverage_RBFN = {}
    negative_d = '/home/leo/Documents/Database/Pipeline_New/Negative_Cores/Sample_0'
    positive_d = '/home/leo/Documents/Database/Pipeline_New/Cores'
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
            all_training_observed_values = Observed_values(training_set, binary)
            
            print('Working on: ' +str(i)+'_'+str(j)+'_0_0_1_2_1perchain')
            all_training_distance_matrix = Distance_matrix(training_set, training_set, square = True)
            
#            all_training_design_matrix = Design_matrix(all_training_distance_matrix)
            
            cross_number = 5            
            percentages = [1, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005]
            for reg in [0]:
    #            percentages = [0.005]
                # Load up the parameter
                parameter = {}
                parameter['all_training_distance_matrix '] = all_training_distance_matrix 
    #            parameter['all_training_design_matrix'] = all_training_design_matrix
                parameter['all_training_observed_values'] = all_training_observed_values 
                parameter['positive_training_set'] = positive_training_set 
                parameter['negative_training_set'] = negative_training_set 
                parameter['training_set'] = training_set 
                parameter['cross_number'] = cross_number 
                parameter['percentages'] = percentages 
                parameter['reg'] = reg
                
                # Cross_validate
                Cross_validation(parameter, method)               
                # Load the results to cross_coverage_RBFN
                key = str(i)+'_'+str(j)+'_0_0_1_2_1perchain_'+str(reg)
                cross_coverage_RBFN[key] = parameter['list_AUC']
        
    return cross_coverage_RBFN
     
######################################################################################
    
def Batch_cross_validation_the_latest(binary, method):
    # The final result will be stored in cross_coverage_RBFN
    cross_coverage_RBFN = {}
    negative_d = '/home/leo/Documents/Database/Pipeline_New/Latest/Negative_cores/sample_0'
    positive_d = '/home/leo/Documents/Database/Pipeline_New/Cores'
    for i in range(1, 5):
        for j in range(1,5):
            p_name_train = 'training_'+str(i)+'_'+str(j)+'_0_0_1_2_1perchain'
            p_name_test = 'testing_'+str(i)+'_'+str(j)+'_0_0_1_2_1perchain'
            os.chdir(positive_d)
            with open(p_name_train, 'r') as f:
                positive_training_set = json.load(f)
            with open(p_name_test, 'r') as f:
                positive_testing_set = json.load(f)
            n_name = p_name_train+'_negative'
            os.chdir(negative_d)
            with open(n_name, 'r') as f:
                negative_training_set = json.load(f)
            
            '''*******The following block is to calculate the observed values 
            ***********according whether the observed should be binary or not'''
            # Calculate 
            training_set = copy.deepcopy(positive_training_set)
            training_set.extend(positive_testing_set)
            training_set.extend(negative_training_set)   
            # Get the observed values
            all_training_observed_values = Observed_values(training_set, binary)
            
            print('Working on: ' +str(i)+'_'+str(j)+'_0_0_1_2_1perchain')
            all_training_distance_matrix = Distance_matrix(training_set, training_set, square = True)
            
#            all_training_design_matrix = Design_matrix(all_training_distance_matrix)
            
            cross_number = 5
            
            percentages = [1, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005]
            for reg in [0]:
                # Load up the parameter
                parameter = {}
                parameter['all_training_distance_matrix '] = all_training_distance_matrix 
    #            parameter['all_training_design_matrix'] = all_training_design_matrix
                parameter['all_training_observed_values'] = all_training_observed_values 
                parameter['positive_training_set'] = positive_training_set 
                parameter['negative_training_set'] = negative_training_set 
                parameter['training_set'] = training_set 
                parameter['cross_number'] = cross_number 
                parameter['percentages'] = percentages 
                parameter['reg'] = reg
                # Cross_validate
                Cross_validation(parameter, method)
                
                
                # Load the results to cross_coverage_RBFN
                key = str(i)+'_'+str(j)+'_0_0_1_2_1perchain_'+str(reg)
                cross_coverage_RBFN[key] = parameter['list_AUC']
        
    return cross_coverage_RBFN
################################################################################
'''
Run the main and Save the results
'''
class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)
###################################################################################
def Cross_validation_in_one_function(method):
    cross_coverage_RBFN_binary = Batch_cross_validation(binary = True, method=method)
    os.chdir('/home/leo/Documents/Database/Pipeline_New/Results')
    with open('cross_binary', 'w') as f:
        json.dump(cross_coverage_RBFN_binary, f, cls=NumpyEncoder)
        
    cross_coverage_RBFN_numerical = Batch_cross_validation(binary = False, method=method)
    os.chdir('/home/leo/Documents/Database/Pipeline_New/Results')
    with open('cross_numerical', 'w') as f:
        json.dump(cross_coverage_RBFN_numerical, f, cls=NumpyEncoder)
    
    cross_coverage_RBFN_the_latest_binary = Batch_cross_validation_the_latest(binary=True, method=method)
    os.chdir('/home/leo/Documents/Database/Pipeline_New/Results')
    with open('cross_latest_binary', 'w') as f:
        json.dump(cross_coverage_RBFN_the_latest_binary, f, cls=NumpyEncoder)
        
    cross_coverage_RBFN_the_latest_numerical = Batch_cross_validation_the_latest(binary=False, method=method)
    os.chdir('/home/leo/Documents/Database/Pipeline_New/Results')
    with open('cross_latest_numerical', 'w') as f:
        json.dump(cross_coverage_RBFN_the_latest_numerical, f, cls=NumpyEncoder)
'''***************************************************************************'''
'''***************************************************************************'''
'''***************************************************************************'''
'''***************************************************************************'''
'''Lets run'''
Cross_validation_in_one_function(method='BFGS')
#########################################################################
'''
From the cross_coverage_RBFN, we see that more the number of the centers, the better
the prediction. Thus we use all the non_redundent_training set as the centers

Test_coverage_RBFN:
    This function is to calculate the precision recall rate for the positive_testing_set
    and the negative_testing_set. There are ten different sets of negative_testing_set. We
    can calculate the separately
Input:
    best_percentage:
        The percentage which performs the best in the cross_coverage_RBFN
Output:
    test_coverage_RBFN_results:
        a list of dictionaries, each dictionary contains the following values
        coeff:
            a list of coefficients
        non_redundant_training_set:
            a sub set of training_set with the redundant removed
        area_list:
            a list of areas, for each negative_testing_set
        recall_list:
            a list of recalls, each recall is a list, corresponding to one 
            set of negative_testing_set.
        precision_list:
            a list of precisions, each precision is a list, corresponding to one
            set of negative_testing_set.
'''
def Test_AUC(binary = True):
    test_coverage_RBFN_results={}
    # Read the training_data and the testing data
    negative_d = '/home/leo/Documents/Database/Pipeline_New/Negative_Cores/Sample_0'
    positive_d = '/home/leo/Documents/Database/Pipeline_New/Cores'
    n_directory = '/home/leo/Documents/Database/Pipeline_New/Negative_Cores/Sample_'
    for i in range(1, 5):
        for j in range(1, 5):
            # set the empty container
            key = str(i)+'_'+str(j)+'_0_0_1_2_1perchain'
            test_coverage_RBFN_results[key] = {}
            
            p_train = 'training_'+str(i)+'_'+str(j)+'_0_0_1_2_1perchain'
            p_test = 'testing_'+str(i)+'_'+str(j)+'_0_0_1_2_1perchain'
            os.chdir(positive_d)
            with open(p_train, 'r') as f:
                positive_training_set = json.load(f)
            with open(p_test, 'r') as f:
                positive_testing_set = json.load(f)
                               
            n_traing = p_train+'_negative'
            n_test = p_test + '_negative'
            os.chdir(negative_d)
            with open(n_traing, 'r') as f:
                negative_training_set = json.load(f)
            print('Working on '+ key)    
            AUC_list, FPR_list, TPR_list, non_redundent_training_set, coeff = \
            Sub_test_AUC(positive_training_set, negative_training_set,positive_testing_set,n_test,n_directory, binary)
                
            test_coverage_RBFN_results[key]['non_redundent_training_set'] = \
                                                    non_redundent_training_set
            test_coverage_RBFN_results[key]['coeff'] = coeff
            test_coverage_RBFN_results[key]['AUC_list'] = AUC_list
            test_coverage_RBFN_results[key]['FPR_list'] = FPR_list
            test_coverage_RBFN_results[key]['TPR_list'] = TPR_list
            
    return test_coverage_RBFN_results
########################################################################
def Sub_test_AUC(positive_training_set, negative_training_set,positive_testing_set,n_test,n_directory, binary):
    # Train the model
    training_set = copy.deepcopy(positive_training_set)
    training_set.extend(negative_training_set)
    # Get the observed values
    all_training_observed_values = Observed_values(training_set, binary)
    all_training_observed_values = np.reshape(all_training_observed_values, (-1,1))

    all_training_distance_matrix = Distance_matrix(training_set, training_set, square = True)            
    all_training_design_matrix = Design_matrix(all_training_distance_matrix)
    
    # Remove the duplicates
    centers, non_redundent_training_set = Remove_duplicates(training_set)
    
    # Extract the design matrix
    centers.append(-1)
    design_matrix = all_training_design_matrix[:,:][:, centers]
    centers.remove(-1)
    
    # Train the model
    # Pay attention to the termination condition
#            print(all_training_observed_values)
    parameter, loss = Train_RBFN_BFGS(design_matrix, all_training_observed_values,\
                                      rho=0.9, c = 1e-3, termination = 1e-2,\
                                     parameter_inheritance = False) 
    # Take out the coeff
    coeff = parameter['coeff']

    
    # Pedict on different testing set
    AUC_list = []
    FPR_list = []
    TPR_list = []
    for n in range(10):
        negative_test_d = n_directory + str(n)
        os.chdir(negative_test_d)
        with open(n_test, 'r') as f:
            negative_testing_set = json.load(f)
        # Reverse the negative testing set
        
        testing_set = copy.deepcopy(positive_testing_set)
        testing_set.extend(negative_testing_set)
        
        # Calculate the observed_values
        observed_values = []
        for test in testing_set:
            if test[2] > 0:
                observed_values.append(1)
            else:
                observed_values.append(test[2])
                
        # Calculate the testing_design_matrix
        print('Length of the testing set:  ', len(testing_set))
        print('Length of the observed values:  ' , len(observed_values))
        print('Calculating the testing design matrix: ', n)
        testing_distance_matrix = Distance_matrix(testing_set, non_redundent_training_set, square = False)
        testing_design_matrix = Design_matrix(testing_distance_matrix)
        
        # Do the prediction
        pred = testing_design_matrix.dot(coeff)
        
        # Calculate the AUC
        AUC, TPR, FPR = Calculate_AUC(pred, observed_values)
#                
#                # Load
        AUC_list.append(AUC)
        FPR_list.append(FPR)
        TPR_list.append(TPR)
    
    return AUC_list, FPR_list, TPR_list, non_redundent_training_set, coeff
#test_coverage_RBFN_results_new = Test_AUC(binary = True)
#test_coverage_RBFN_results_new['4_4_0_0_1_2_1perchain']['coeff'][:6]
#test_coverage_RBFN_results['4_4_0_0_1_2_1perchain']['coeff'][:6]
############################################################################
def Test_AUC_the_latest(binary = True):
    test_coverage_RBFN_results={}
    # Read the training_data and the testing data
    negative_d = '/home/leo/Documents/Database/Pipeline_New/Latest/Negative_cores/sample_0'
    positive_d = '/home/leo/Documents/Database/Pipeline_New/Cores'
    positive_d_latest = '/home/leo/Documents/Database/Pipeline_New/Latest/cores'
    n_directory = '/home/leo/Documents/Database/Pipeline_New/Latest/Negative_cores/sample_'
    for i in range(1, 5):
        for j in range(1, 5):
            # set the empty container
            key = str(i)+'_'+str(j)+'_0_0_1_2_1perchain'
            test_coverage_RBFN_results[key] = {}
            
            p_train = 'training_'+str(i)+'_'+str(j)+'_0_0_1_2_1perchain'
            p_test = 'testing_'+str(i)+'_'+str(j)+'_0_0_1_2_1perchain'
            os.chdir(positive_d)
            with open(p_train, 'r') as f:
                positive_training = json.load(f)
#            with open(p_test, 'r') as f:
#                positive_testing = json.load(f)
            positive_training_set = copy.deepcopy(positive_training)
#            positive_training_set.extend(positive_testing)
                               
            n_traing = p_train+'_negative'
            n_test = p_test + '_negative'
            os.chdir(negative_d)
            with open(n_traing, 'r') as f:
                negative_training_set = json.load(f)
                
            os.chdir(positive_d_latest)
            with open(p_test, 'r') as f:
                positive_testing_set = json.load(f)
                
            print('Working on '+ key)    
            AUC_list, FPR_list, TPR_list, non_redundent_training_set, coeff = \
            Sub_test_AUC(positive_training_set, negative_training_set,positive_testing_set,n_test,n_directory, binary)
                
            test_coverage_RBFN_results[key]['non_redundent_training_set'] = \
                                                    non_redundent_training_set
            test_coverage_RBFN_results[key]['coeff'] = coeff
            test_coverage_RBFN_results[key]['AUC_list'] = AUC_list
            test_coverage_RBFN_results[key]['FPR_list'] = FPR_list
            test_coverage_RBFN_results[key]['TPR_list'] = TPR_list
            
    return test_coverage_RBFN_results

'''
Run the above function
'''
#test_coverage_RBFN_results = Calculate_recall_precision(binary = True)

#######################################################################
def Test_AUC_in_one_function():
    test_coverage_RBFN_results_binary = Test_AUC(binary = True)
    os.chdir('/home/leo/Documents/Database/Pipeline_New/Results')
    with open('test_coverage_RBFN_results_binary', 'w') as f:
        json.dump(test_coverage_RBFN_results_binary , f, cls=NumpyEncoder)
        
    test_coverage_RBFN_results_numerical = Test_AUC(binary = False)
    os.chdir('/home/leo/Documents/Database/Pipeline_New/Results')
    with open('test_coverage_RBFN_results_numerical', 'w') as f:
        json.dump(test_coverage_RBFN_results_numerical, f, cls=NumpyEncoder)
        
    test_coverage_RBFN_results_latest_binary = Test_AUC_the_latest(binary = True)
    os.chdir('/home/leo/Documents/Database/Pipeline_New/Results')
    with open('test_coverage_RBFN_results_latest_binary', 'w') as f:
        json.dump(test_coverage_RBFN_results_latest_binary, f, cls=NumpyEncoder)
        
    test_coverage_RBFN_results_latest_numerical = Test_AUC_the_latest(binary = False)
    os.chdir('/home/leo/Documents/Database/Pipeline_New/Results')
    with open('test_coverage_RBFN_results_latest_numerical', 'w') as f:
        json.dump(test_coverage_RBFN_results_latest_numerical, f, cls=NumpyEncoder)
    # Save the results




    
#####################################################################

'''
Reverse_test:
    A function to test the testing set with the Ab and Ag sequences switched
Input:
    test_coverage_RBFN_results
Output:
    reverse_test_coverage_RBFN_results
    ['recalls_list']
    [precidions_list]
    [area_list]
'''
def Discrimation_test(test_coverage_RBFN_results):
    reverse_test_coverage_RBFN_results={}
    # Read the training_data and the testing data
    positive_d = '/home/leo/Documents/Database/Pipeline_New/Cores'
    for i in range(1, 5):
        for j in range(1, 5):
            # set the empty container
            key = str(i)+'_'+str(j)+'_0_0_1_2_1perchain'
            reverse_test_coverage_RBFN_results[key] = {}
            non_redundent_training_set = test_coverage_RBFN_results[key]['non_redundent_training_set']
            coeff = test_coverage_RBFN_results[key]['coeff']            

            p_test = 'testing_'+str(i)+'_'+str(j)+'_0_0_1_2_1perchain'
            os.chdir(positive_d)
            with open(p_test, 'r') as f:
                positive_testing_set = json.load(f)            
            
            n_test = 'testing_'+str(j)+'_'+str(i)+'_0_0_1_2_1perchain'
            with open(n_test, 'r') as f:
                negative_testing = json.load(f)
#            n_train = 'training_'+str(j)+'_'+str(i)+'_0_0_1_2_1perchain'
#            with open(n_train, 'r') as f:
#                negative_training = json.load(f)
#            negative_testing.extend(negative_training)
            # Reverse the negative testing set
            negative_testing_set = []
            for parepi_negative in negative_testing:
                negative_testing_set.append([parepi_negative[1], parepi_negative[0],\
                                             -1, parepi_negative[3]])                    
            
            testing_set = copy.deepcopy(positive_testing_set)
            testing_set.extend(negative_testing_set)
            
            # Calculate the observed_values
            observed_values = []
            for test in testing_set:
                if test[2] > 0:
                    observed_values.append(1)
                else:
                    observed_values.append(test[2])
                    
            # Calculate the testing_design_matrix
            print('Length of the testing set:  ', len(testing_set))
            print('Length of the observed values:  ' , len(observed_values))
            print('Calculating the testing design matrix: ')
            testing_distance_matrix = Distance_matrix(testing_set, non_redundent_training_set, square = False)
            testing_design_matrix = Design_matrix(testing_distance_matrix)
            
            # Do the prediction
            pred = testing_design_matrix.dot(coeff)

            AUC, TPR, FPR = Calculate_AUC(pred, observed_values)
            
                
            reverse_test_coverage_RBFN_results[key]['AUC'] = AUC
            reverse_test_coverage_RBFN_results[key]['TPR'] = TPR
            reverse_test_coverage_RBFN_results[key]['FPR'] = FPR
            
    return reverse_test_coverage_RBFN_results

#######################################################################
def Discrimation_test_the_latest(test_coverage_RBFN_results):
    reverse_test_coverage_RBFN_results={}
    # Read the training_data and the testing data
    positive_d = '/home/leo/Documents/Database/Pipeline_New/Latest/cores'
    for i in range(1, 5):
        for j in range(1, 5):
            # set the empty container
            key = str(i)+'_'+str(j)+'_0_0_1_2_1perchain'
            reverse_test_coverage_RBFN_results[key] = {}
            non_redundent_training_set = test_coverage_RBFN_results[key]['non_redundent_training_set']
            coeff = test_coverage_RBFN_results[key]['coeff']            

            p_test = 'testing_'+str(i)+'_'+str(j)+'_0_0_1_2_1perchain'
            os.chdir(positive_d)
            with open(p_test, 'r') as f:
                positive_testing_set = json.load(f)            
            # Let the training set of the match-type(j, i) be the negative training set
            n_test = 'testing_'+str(j)+'_'+str(i)+'_0_0_1_2_1perchain'
            n_train = 'training_'+str(j)+'_'+str(i)+'_0_0_1_2_1perchain'
            with open(n_test, 'r') as f:
                negative_testing = json.load(f)
            os.chdir('/home/leo/Documents/Database/Pipeline_New/Cores')
            with open(n_test, 'r') as f:
                negative_testing_old = json.load(f)
            with open(n_train, 'r') as f:
                negative_training_old = json.load(f)
            # Reverse the negative testing set
            # Combine all of them
            negative_testing.extend(negative_testing_old)
            negative_testing.extend(negative_training_old)
            negative_testing_set = []
            for parepi_negative in negative_testing:
                negative_testing_set.append([parepi_negative[1], parepi_negative[0],\
                                             -1, parepi_negative[3]])                    
            
            testing_set = copy.deepcopy(positive_testing_set)
            testing_set.extend(negative_testing_set)
            
            # Calculate the observed_values
            observed_values = []
            for test in testing_set:
                if test[2] > 0:
                    observed_values.append(1)
                else:
                    observed_values.append(test[2])
                    
            # Calculate the testing_design_matrix
            print('Length of the testing set:  ', len(testing_set))
            print('Length of the observed values:  ' , len(observed_values))
            print('Calculating the testing design matrix: ')
            testing_distance_matrix = Distance_matrix(testing_set, non_redundent_training_set, square = False)
            testing_design_matrix = Design_matrix(testing_distance_matrix)
            
            # Do the prediction
            pred = testing_design_matrix.dot(coeff)

            AUC, TPR, FPR = Calculate_AUC(pred, observed_values)
            
                
            reverse_test_coverage_RBFN_results[key]['AUC'] = AUC
            reverse_test_coverage_RBFN_results[key]['TPR'] = TPR
            reverse_test_coverage_RBFN_results[key]['FPR'] = FPR
            
    return reverse_test_coverage_RBFN_results
          
#########################################################################3
def Reverse_test_in_one_function():
#    os.chdir('/home/leo/Documents/Database/Pipeline_New/Codes/Results')
#    with open('test_coverage_RBFN_results_binary', 'r') as f:
#        test_coverage_RBFN_results_binary = json.load(f)
#    reverse_test_binary = \
#                    Discrimation_test(test_coverage_RBFN_results_binary)
    
    os.chdir('/home/leo/Documents/Database/Pipeline_New/Codes/Results')
    with open('test_coverage_RBFN_results_latest_binary', 'r') as f:
        test_coverage_RBFN_results_latest_binary = json.load(f)
    reverse_test_latest_binary = \
                    Discrimation_test_the_latest(test_coverage_RBFN_results_latest_binary)
                    
#    return reverse_test_binary, reverse_test_latest_binary
    return reverse_test_latest_binary
        
##########################################################################################
def Selected_centers(cut = 6):
    os.chdir('/home/leo/Documents/Database/Pipeline_New/Codes/Results')
    with open('one_to_one_frequency', 'r') as f:
        one_to_one_frequency = json.load(f)

    selected_positive_cores = []
    for one in one_to_one_frequency:
        if one[2]>= 6:
            selected_positive_cores.append([[one[0]],[ one[1]], 1])
    return selected_positive_cores
#selected_positive_cores = Selected_centers(cut = 6)
#with open('selected_one_to_one', 'w') as f:
#    json.dump(selected_positive_cores, f)
###############################################################################
def Test_AUC_selected_centers(binary = True):
    test_coverage_RBFN_results={}
    # Read the training_data and the testing data
    negative_d = '/home/leo/Documents/Database/Pipeline_New/Negative_Cores/Sample_0'
    positive_d = '/home/leo/Documents/Database/Pipeline_New/Cores'
    n_directory = '/home/leo/Documents/Database/Pipeline_New/Negative_Cores/Sample_'
    for i in range(1, 2):
        for j in range(1, 2):
            # set the empty container
            key = str(i)+'_'+str(j)+'_0_0_1_2_1perchain'
            test_coverage_RBFN_results[key] = {}
            
            p_test = 'testing_'+str(i)+'_'+str(j)+'_0_0_1_2_1perchain'
            os.chdir(positive_d)
            with open('selected_one_to_one', 'r') as f:
                positive_training_set = json.load(f)
            with open(p_test, 'r') as f:
                positive_testing_set = json.load(f)
                               
            n_traing ='training_'+ key+'_negative'
            n_test = p_test + '_negative'
            os.chdir(negative_d)
            with open(n_traing, 'r') as f:
                negative_training_set = json.load(f)
            print('Working on '+ key)    
            AUC_list, FPR_list, TPR_list, non_redundent_training_set, coeff = \
            Sub_test_AUC(positive_training_set, negative_training_set,positive_testing_set,n_test,n_directory, binary)
                
            test_coverage_RBFN_results[key]['non_redundent_training_set'] = \
                                                    non_redundent_training_set
            test_coverage_RBFN_results[key]['coeff'] = coeff
            test_coverage_RBFN_results[key]['AUC_list'] = AUC_list
            test_coverage_RBFN_results[key]['FPR_list'] = FPR_list
            test_coverage_RBFN_results[key]['TPR_list'] = TPR_list
            
    return test_coverage_RBFN_results
##################################################################################
#reverse_test_latest_binary_all_reverse_negative = Reverse_test_in_one_function()
#reverse_test_latest_binary_all_reverse_negative['4_1_0_0_1_2_1perchain'].keys()

#os.chdir('/home/leo/Documents/Database/Pipeline_New/Codes/Results')
#with open('reverse_test_latest_binary_all_negative', 'w') as f:
#    json.dump(reverse_test_latest_binary_all_reverse_negative, f)

#selected_test_results = Test_AUC_selected_centers(binary = True)
#os.chdir('/home/leo/Documents/Database/Pipeline_New/Codes/Results')
#with open('selected_test_results', 'w') as f:
#    json.dump(selected_test_results, f, cls=NumpyEncoder)
#if __name__ == '__main__':
#    Cross_validation_in_one_function()
#Test_AUC_in_one_function()
#reverse_test_binary, reverse_test_latest_binary = Reverse_test_in_one_function()
#os.chdir('/home/leo/Documents/Database/Pipeline_New/Results')
#with open('test_coverage_RBFN_results_binary', 'r') as f:
#    test_coverage_RBFN_results_binary = json.load(f)
#os.chdir('/home/leo/Documents/Database/Pipeline_New/Results')
#with open('reverse_test_binary', 'w') as f:
#    json.dump(reverse_test_binary, f)
#with open('reverse_test_latest_binary', 'w') as f:
#    json.dump(reverse_test_latest_binary, f)

#reverse_test_binary['2_1_0_0_1_2_1perchain']['AUC']
#FPR = reverse_test_binary['2_1_0_0_1_2_1perchain']['FPR']
#TPR = reverse_test_binary['2_1_0_0_1_2_1perchain']['TPR']
#plt.plot(FPR, TPR)
#plt.plot([0,1], [0,1])
#with open('test_coverage_RBFN_results_latest_binary', 'r') as f:
#    test_coverage_RBFN_results_latest_binary =  json.load(f)
#with open('test_coverage_RBFN_results_numerical', 'r') as f:
#    test_coverage_RBFN_results_numerical = json.load(f)
#with open('test_coverage_RBFN_results_latest_numerical', 'r') as f:
#    test_coverage_RBFN_results_latest_numerical =  json.load(f)
##    
#test_coverage_RBFN_results_latest_binary['2_2_0_0_1_2_1perchain'].keys()
#FPR_list = test_coverage_RBFN_results_latest_binary['2_2_0_0_1_2_1perchain']['FPR_list']
#TPR_list = test_coverage_RBFN_results_latest_binary['2_2_0_0_1_2_1perchain']['TPR_list']
#for i in range(len(FPR_list)):
#    FPR = FPR_list[i]
#    TPR =TPR_list[i]
#    plt.plot(FPR, TPR)
#plt.plot([0,1], [0,1])
#np.average(np.array(test_coverage_RBFN_results_latest_binary['2_2_0_0_1_2_1perchain']['AUC_list']))
#test_coverage_RBFN_results_latest_binary['2_2_0_0_1_2_1perchain']['AUC_list']
#test_coverage_RBFN_results_numerical['2_2_0_0_1_2_1perchain']['AUC_list']
#test_coverage_RBFN_results_latest_numerical['2_2_0_0_1_2_1perchain']['AUC_list']



