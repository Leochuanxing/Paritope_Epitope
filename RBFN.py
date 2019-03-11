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
os.chdir('/home/leo/Documents/Database/Pipeline/Codes and testing data')
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
####################################################################
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
    Ab_seq1 = parepi_1[0]
    Ab_seq2 = parepi_2[0]
    Ag_seq1 = parepi_1[1]
    Ag_seq2 = parepi_2[1]
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
###############################################################################
'''
Loss:
    To return the loss, the gradient
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
############################################################################
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
#############################################################################
'''
The purpose of this block is to do cross validation, select centers and make predictiion about the 
testing data accroding to the results of the cross validation

parameter contains
    parameter['percentages']
    parameter['training_design_matrix']
    parameter['training_observed_values']
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
    parameter['centers']
    parameter['coeff']
    parameter['reg']
    
******************Explanation about the above values***********************
percentages:
    a list of percentages, which are the hyperparameters to be selected by the cross validation
    
training_design_matrix:
    is the design matrix with the training set as the rows and the same training set as the centers
    
training_observed_values:
    a list of observed values extracted from the training set, and the orders are corresponding with each other
    
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
################################################################################################

'''
One_step_reduce_centers:
    a function to reduce the number of centers by 1
Input:
    parameter:
        a dictionary, contains all the parameters related to the model
        
        parameter['centers']:
            a list of centers, this list will be reduced by 1 element each time
        parameter['design_matrix']:
            a design matrix, with the columns corresponding to the centers
                    
      ##### Pay attention that the column number is one larger than the row number.#####
            
        parameter['observed_values']:
            an array in the shape of (np.shape(design_matrix), ), with length equal to the number of the rows of the design_matrix
        parameter['coeff']:
            An array in the shape of (len(centers)+1, 1), gives the coefficient of 
            corresponding to the centers and the constant term
        parameter['reg']:
            an array in the shape of (len(centers)+1, 1), gives the coefficient of the
            regularization term.

Output:
    parameter, a dictionary with all the values updated
'''
'''
We can add a control to speed up the process
'''
def One_step_reduce_centers(parameter):
    # Take out the values
    centers = parameter['centers']
    design_matrix = parameter['design_matrix']
    training_observed_values = parameter['training_observed_values']

#    ratio = 3000/len(centers)
    if len(centers) > 1500 :
        termination = 1.5*len(centers)
    elif len(centers) <= 1500 and len(centers) >1000:
        termination =  len(centers)
    elif len(centers)<= 1000 and len(centers)>500: 
        termination = 10
    elif len(centers) <= 500:
        termination = len(centers)/1000
#    termination =10* len(centers)

        
    parameter, loss = Train_RBFN_BFGS(design_matrix, training_observed_values, rho=0.85, c = 1e-3, termination=termination,\
                                      parameter_inheritance = True, parameter=parameter)
    
    # Take out the coefficient and get rid of the unimportant centers
    coeff = parameter['coeff']
    # Take the absolute value
    coeff_list = np.abs(coeff)
    # match up the absolute values of coeff with other indices
    match_up = []
    for i in range(len(centers)):
        match_up.append([centers[i], coeff_list[i], i])
    # sort the match up according to the value of coeff
    match_up.sort(key = lambda x:x[1])
    # To speed up the process, we remove more than one centers if the number of the centers is large.
    # Creat an empty container to contain the centers to be removed.
    removed_sets = []
    if len(match_up) > 2000:
        n_remove = math.floor(len(centers)*0.004)
    elif len(match_up) > 1500 and len(match_up) <= 2000:
        n_remove = math.floor(len(centers)*0.002)
    else:
        n_remove = 1
    # Load up the removed_sets according to the number of removed centers
    for i in range(n_remove):
        removed_sets.append(match_up[i])
    
    # Update the centers
    removed_col = []
    for to_be_removed in removed_sets:    
        centers.remove(to_be_removed[0])
        removed_col.append(to_be_removed[2])
    # Update the coefficient
    coeff = np.delete(coeff, removed_col)
    coeff = np.reshape(coeff, (-1, 1))
    # Update the design matrix
    design_matrix = np.delete(design_matrix, removed_col, 1)
    # Load the updated values to the parameter
    parameter['centers']= centers
    parameter['coeff'] = coeff
    parameter['reg'] = np.ones((len(coeff), 1))
    parameter['design_matrix'] = design_matrix
    return parameter

###############################################################################

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

#centers = Remove_duplicates(training_set)

#############################################################

'''
Generate_coeff_centers:
    a function to generate the coeff centers with the length of the centers equal to a given length
Input:
    parameter:
        a dictionary with one more value
        parameter['list_n_centers']: a list of integers, gives the number of centers to be kept
#    control_coeff:
#        a float, to speed up the selection of centers, when removing the centers, we remove the ones with 
#        the absolute values fo the coefficients less than control_coeff, otherwise, we remove the one with 
#        the least absolute value of the coefficient.
Output:
    parameter:
        with one more value
        parameter['list_centers']: a list with each element a list of centers
        parameter['list_coeff']: a list of coeff, each coeff is corresponding to a list of centers
        parameter['list_reg']: similar as above
'''

def Coeff_select_centers(parameter):
    # Take out the values from the parameter
    list_n_centers = parameter['list_n_centers']
    centers = parameter['centers']
    n_centers = len(centers)
    # Create empty lists
    list_centers = []
    list_coeff = []
    list_reg = []
    # Reduce the centers according the list_n_centers
    for n in list_n_centers:
        while n < n_centers:
            parameter = One_step_reduce_centers(parameter)
            # update the n_centers
            centers = parameter['centers']
            n_centers = len(centers)
            print('Length of centers: ', n_centers)
        list_centers.append(copy.deepcopy(centers))
        # Load list_coeff and list_reg
        list_coeff.append(copy.deepcopy(parameter['coeff']))
        list_reg.append(copy.deepcopy(parameter['reg']))
    
    # Load list_centers, list_coeff and list_reg to the parameter
    parameter['list_centers'] = list_centers
    parameter['list_coeff'] = list_coeff
    parameter['list_reg'] = list_reg
    
    return parameter
##################################################################################################
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
Cross_validation:
    a function to caltulate the area under the precision_recall curve, this function is specifically 
    designed for the cross_validation
Input:
    parameter['training_set'], parameter[testing_set]
Output:
    parameter with one more value
    parameter['list_recall_precision'], a list of recall list and precision list
'''
def Cross_validation(parameter):

    # Take out the values
    training_design_matrix = parameter['training_design_matrix']
    training_observed_values = parameter['training_observed_values']# pay attention to this one, it is used in the one step reduce
    positive_training_set = parameter['positive_training_set']
    negative_training_set = parameter['negative_training_set']
    cross_number = parameter['cross_number']
    percentages = parameter['percentages']
    
    # Initiate list_area_under_rp_curve 
    list_area_under_rp_curve = [0*i for i in range(cross_number)]
        
    # Generate the cross indices
    cross_train_indices, cross_test_indices = \
    Raised_cross_indices(positive_training_set, negative_training_set, cross_number)
    
    
    # Go into the cross
    for i_cross in range(cross_number):
        # Load the initial design matrix
        row_indices = cross_test_indices[i_cross]
        col_indices= copy.deepcopy(cross_train_indices[i_cross])
        col_indices.append(-1)# attach the last column of the constant 1
        parameter['design_matrix'] = training_design_matrix[row_indices:][:col_indices]
        # Load the list_n_centers
        parameter['list_n_centers'] = []
        for per in percentages:
            n_center = math.floor(len(cross_train_indices[i_cross])*per)
            if n_center >= 1:
                parameter['list_n_centers'].append(n_center)

                

    # Calculate the big_testing_design_matrix
    big_training_set = []
    for i in big_centers:
        big_training_set.append(training_set [i])
    big_test_design_matrix = Testing_design_matrix(testing_set, big_training_set)
    # Go into one set of centers
    for i in range(len(list_centers)):
        centers_pos = list_relative_pos[i]
        coeff = list_coeff[i]
        centers_pos.append(-1)
        test_design_matrix = big_test_design_matrix[:, centers_pos]
        
        pred = test_design_matrix.dot(coeff)
        # Match up the observed values with the pred
        match_up = []
        for i in range(len(observed_values)):
            match_up.append([pred[i], observed_values[i]])
        match_up.sort(key = lambda x:x[0], reverse = True)
            
        
        # Calculate the denominator
        denominator = 0
        for i in observed_values:
            if i == 1:
                denominator += 1
        # Calculate the recall and precision        
        n_positive = 0
        n_negative = 0
        precision = []
        recall = []
        for match in match_up:
            if match[1] == 1:
                n_positive += 1
            else:
                n_negative += 1
            precision.append(n_positive/(n_positive + n_negative))
            recall.append(n_positive/denominator)
        list_recall_precision.append([copy.deepcopy(recall), copy.deepcopy(precision)])
        
    parameter['list_recall_precision'] = list_recall_precision
    pass
########################################################################################3    

  
    
#####################################################################################
# Record the time

def main(wd_results, wd_negative_samples, header):
    os.chdir(wd_results)
    with open(header+'_aa_train', 'r') as f:
        positive_training_set = json.load(f)
    with open(header+'_aa_test', 'r') as f:
        positive_testing_set = json.load(f)
        
    os.chdir(wd_negative_samples)
    with open(header+'_train_negative', 'r') as f:
        negative_training_set = json.load(f)
    with open(header+'_test_negative', 'r') as f:
        negative_testing_set = json.load(f)  
    
    # set the observed value as 1 for the positive trainging set.    
    for parepi in positive_training_set:
        parepi[2] = 1
    for parepi in positive_testing_set:
        parepi[2] = 1
    

    
    training_set = copy.deepcopy(positive_training_set)
    training_set.extend(negative_training_set)
    # Remove the duplicates
    centers, non_redundent_training_set = Remove_duplicates(training_set)
    # Calculate the distance matrix and the total design matrix
    print('Calculating the distance matrix')
    distance_matrix = Distance_matrix(training_set, mode='Multiplication') # mode is an ajustable parameter
    design_matrix = Design_matrix(distance_matrix, basis_function = 'Markov', radius = 1) # basis function is adjustable

    # Remove the duplicates
    # Load the parameter
    parameter = {}
    parameter['design_matrix'] = design_matrix
    parameter['centers'] = centers
    parameter['non_redundent_training_set'] = non_redundent_training_set
    parameter['training_set'] = training_set
    len(centers)
    # Initiate the coeff and reg
    parameter['coeff'] = np.zeros((len(centers)+1,1))
    parameter['reg'] = np.ones((len(centers)+1,1))
        
    # Generate the observed values
    observed_values = np.zeros((len(training_set ), 1))
    for i in range(len(training_set)):
        observed_values[i,0]=training_set [i][2]
    # Load the parameter
    parameter['observed_values'] = observed_values

    # Generate the list_n_centers
    percentage = [0.5, 0.2,0.1,0.05,0.02,0.01,0.008,0.005,0.002,0.001]
    list_n_centers = []
    for i in percentage:
        n= math.floor(len(non_redundent_training_set)*i)
        if n>= 1:
            list_n_centers.append(n)
        
    parameter['list_n_centers'] = list_n_centers
    
    parameter = Coeff_select_centers(parameter)



    testing_set = copy.deepcopy(positive_testing_set)
    testing_set.extend(negative_testing_set)
    observed_values = []
    for i in testing_set:
        observed_values.append(i[2])
    Recall_precision_calculation(parameter, testing_set, observed_values)
    
    return parameter
#len(parameter['list_recall_precision'])
#list_recall_precision = parameter['list_recall_precision']
#os.chdir("/home/leo/Documents/Database/Pipeline/Results/1_free")
#with open('2_2_test_results')as f:
#    coverage_results = json.load(f)
#########################################################################
class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)
######################################################
ordered_headers = ['2_2','3_3', '2_3', '3_2', '3_4', '4_3','2_4', '4_2','4_4',\
                   '1_1','1_2', '2_1', '1_3', '3_1', '1_4', '4_1']
############################################################
if __name__ == '__main__':
    directory_paires = [['/home/leo/Documents/Database/Pipeline/Results/1_free',\
      '/home/leo/Documents/Database/Pipeline/Negative samples/1_free'],\
     ['/home/leo/Documents/Database/Pipeline/Results/0_free',\
                           '/home/leo/Documents/Database/Pipeline/Negative samples/0_free']]
    
    # Let check if all the corresponding length of the positive training and the negative training
    # are the same
    for i in range(2):
        wd_results = directory_paires[i][0]
        wd_negative_samples = directory_paires[i][1]
        n = 0
        for i in range(1,5):
            for j in range(1, 5):
                header = str(i)+'_'+str(j)
                os.chdir(wd_results)
                with open(header+'_aa_train', 'r') as f:
                    positive_training_set = json.load(f)
                    
                os.chdir(wd_negative_samples)
                with open(header+'_train_negative', 'r') as f:
                    negative_training_set = json.load(f)
                
                if len(positive_training_set) != len(negative_training_set):
                    print(header+' positive and negative are not of the same length')
                    print(wd_results)
                    n = 1
                    break
    if n == 0:
        print('Now, we are going to the large amount of data processing. It may '+\
              'spend a lot of time.')
            
    # Get into the work
    # Fist let us set the working oder of different files. working on the most interesting file 
    # first then gradually reduce to the most uninteresting file, so the even if there are some
    # problems pop up in the middle, we have finished the most interesting files.
    for i in range(2):
        wd_results = directory_paires[i][0]
        wd_negative_samples = directory_paires[i][1]
        for header in ordered_headers:
            print('Working on '+ header)
            results =  main(wd_results, wd_negative_samples, header)
            
             # Save the results
            os.chdir(wd_results)
            with open(header+'_test_results_Coeff_RBFN', 'w') as f:
                json.dump(results, f,  cls=NumpyEncoder)
############################################################
#os.chdir("/home/leo/Documents/Database/Pipeline/Results/1_free")
#with open('2_2_test_results_coeff_RBFN','r') as f:
#    test_results_coeff_RBFN = json.load(f)
#test_results_coeff_RBFN.keys()
#coverage_recall = test_results_coeff_RBFN['recall']
#list_recall_precision = test_results_coeff_RBFN['list_recall_precision']
#len(test_results_coeff_RBFN['centers'])
#for i in test_results_coeff_RBFN['list_n_centers']:
#    print(i)
#for i in range(len(list_recall_precision)):
#    plt.plot(list_recall_precision[i][0], list_recall_precision[i][1])
#    plt.show()
#####################################################################################





