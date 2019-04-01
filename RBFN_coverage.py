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
                if match[1] > 0:
                    '''Here we use >0 to make sure it works for both binary and non binary'''
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

################################################################################
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
####################################################################
def Observed_values(training_set, binary):
    all_training_observed_values = []
    if binary:
        for train in training_set:
            if train[2] >0:
                all_training_observed_values.append(1)
            else:
                all_training_observed_values.append(train[2])
        
        all_training_observed_values = np.array(all_training_observed_values)
    
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
                all_training_observed_values.append(train[2]/rang)
            else:
                all_training_observed_values.append(train[2])
        
        all_training_observed_values = np.array(all_training_observed_values)
        
        return all_training_observed_values
    
###################################################################################
def Batch_cross_validation(binary = True):
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
            key = str(i)+'_'+str(j)+'_0_0_1_2_1perchain'
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
#cross_coverage_RBFN_0_0 = Batch_cross_validation()

# Save the results
#os.chdir('/home/leo/Documents/Database/Pipeline_New/Results')
#with open('cross_coverage_RBFN_0_0', 'w') as f:
#    json.dump(cross_coverage_RBFN_0_0, f)

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
def Calculate_recall_precision(binary = True):
    test_coverage_RBFN_results={}
    # Read the training_data and the testing data
    negative_d = '/home/leo/Documents/Database/Pipeline_New/Negative_Cores/Sample_0'
    positive_d = '/home/leo/Documents/Database/Pipeline_New/Cores'
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
                
            # Train the model
            training_set = copy.deepcopy(positive_training_set)
            training_set.extend(negative_training_set)
            # Get the observed values
            all_training_observed_values = Observed_values(training_set, binary)
            all_training_observed_values = np.reshape(all_training_observed_values, (-1,1))

            print('Calculating the design matrix')
            all_training_distance_matrix = Distance_matrix(training_set, training_set, square = True)            
            all_training_design_matrix = Design_matrix(all_training_distance_matrix)
            
            # Remove the duplicates
            centers, non_redundent_training_set = Remove_duplicates(training_set)
            # Load
            test_coverage_RBFN_results[key]['non_redundent_training_set'] = non_redundent_training_set
            
            # Extract the design matrix
            centers.append(-1)
            design_matrix = all_training_design_matrix[:,:][:, centers]
            centers.remove(-1)
            
            # Train the model
            # Pay attention to the termination condition
            parameter, loss = Train_RBFN_BFGS(design_matrix, all_training_observed_values,\
                                              rho=0.9, c = 1e-3, termination = 1e-2,\
                                             parameter_inheritance = False) 
            # Take out the coeff
            coeff = parameter['coeff']
            # Load
            test_coverage_RBFN_results[key]['coeff'] = coeff
            
            # Pedict on different testing set
            area_list = []
            precisions_list = []
            recalls_list = []
            for n in range(10):
                negative_test_d = '/home/leo/Documents/Database/Pipeline_New/Negative_Cores/Sample_' + str(n)
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

                # Match up the observed values with the pred
                match_up = []
                for j in range(len(observed_values)):
                    match_up.append([pred[j], observed_values[j]])
                match_up.sort(key = lambda x:x[0], reverse = True)
                
                # Calculate the area, precision, recall and precision
                area = 0
                n_positive = 0
                n_negative = 0
                denominator = len(positive_testing_set)
                recalls = []
                precisions = []
                for match in match_up:
                    if match[1] > 0:
                        n_positive += 1
                    elif match[1] == -1:
                        n_negative += 1
                    precision = n_positive/(n_positive+n_negative)
                    area += precision
                    recall = n_positive/denominator
                    
                    recalls.append(recall)
                    precisions.append(precision)                    
                
                # Take the average
                area/=(len(testing_set))
                
                # Load
                area_list.append(area)
                recalls_list.append(recalls)
                precisions_list.append(precisions)
                
            test_coverage_RBFN_results[key]['area_list'] = area_list
            test_coverage_RBFN_results[key]['recalls_list'] = recalls_list
            test_coverage_RBFN_results[key]['precisions_list'] = precisions_list
            
    return test_coverage_RBFN_results
#######################################################################
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
def Reverse_test(test_coverage_RBFN_results):
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
            

            p_test = 'testing_'+str(j)+'_'+str(i)+'_0_0_1_2_1perchain'
            os.chdir(positive_d)
            with open(p_test, 'r') as f:
                positive_testing = json.load(f)            
            # Reverse the testing set
            positive_testing_set = []
            for parepi in positive_testing:
                positive_testing_set.append([parepi[1], parepi[0], parepi[2],parepi[3]])
                            
            # Pedict on different testing set
            area_list = []
            precisions_list = []
            recalls_list = []
            for n in range(10):
                negative_test_d = '/home/leo/Documents/Database/Pipeline_New/Negative_Cores/Sample_' + str(n)
                os.chdir(negative_test_d)
                n_test = p_test + '_negative'
                with open(n_test, 'r') as f:
                    negative_testing = json.load(f)
                # Reverse the negative testing set
                negative_testing_set = []
                for parepi_negative in negative_testing:
                    negative_testing_set.append([parepi_negative[1], parepi_negative[0],\
                                                 parepi_negative[2], parepi_negative[3]])                    
                
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

                # Match up the observed values with the pred
                match_up = []
                for j in range(len(observed_values)):
                    match_up.append([pred[j], observed_values[j]])
                match_up.sort(key = lambda x:x[0], reverse = True)
                
                # Calculate the area, recall and precision
                area = 0
                n_positive = 0
                n_negative = 0
                denominator = len(positive_testing_set)
                recalls = []
                precisions = []
                for match in match_up:
                    if match[1] > 0:
                        n_positive += 1
                    elif match[1] == -1:
                        n_negative += 1
                    precision = n_positive/(n_positive+n_negative)
                    area += precision
                    recall = n_positive/denominator
                    
                    recalls.append(recall)
                    precisions.append(precision)                    
                
                # Take the average
                area/=(len(testing_set))
                
                # Load
                area_list.append(area)
                recalls_list.append(recalls)
                precisions_list.append(precisions)
                
            reverse_test_coverage_RBFN_results[key]['area_list'] = area_list
            reverse_test_coverage_RBFN_results[key]['recalls_list'] = recalls_list
            reverse_test_coverage_RBFN_results[key]['precisions_list'] = precisions_list
            
    return reverse_test_coverage_RBFN_results
    
#######################################################################
          
              
#    os.chdir('/home/leo/Documents/Database/Pipeline_New/Results')    
#        
#    with open('test_coverage_RBFN_0_0_results', 'w') as f:
#        json.dump(test_coverage_RBFN_0_0_results, f, cls = NumpyEncoder)
#    with open('reverse_test_coverage_RBFN_0_0_results', 'w') as f:
#        json.dump(reverse_test_coverage_RBFN_0_0_results , f, cls = NumpyEncoder)
 


#with open('reverse_test_coverage_RBFN_0_0_results', 'r') as f:
#    res_00 = json.load(f)
#    
#average1=[]
#for key, value in res_00.items():
#    average1.append(np.average(np.array(value['area_list'])))
##
#with open('test_coverage_RBFN_results', 'r') as f:
#    res = json.load(f)
#    
#average2=[]    
#for key, value in res.items():
#    average2.append(np.average(np.array(value['area_list'])))
#np.average(np.array(average1)-np.array(average2))

#test_coverage_RBFN_results['1_1_1_1_1_2_1perchain']['recalls_list'][0][-50:]
##########################################################################
'''
Recall_to_ROC_AUC:
    This fucntion is to calculate the ROC and the AUC accroding to the recall list
Input:
    recall:
        a list of recall rate
Output:
    FPR: a list of false positive rate
    TPR: a list of true positive rate
    AUC: a float, gives the area under the curve
'''
def Recall_to_ROC_AUC(recall):
    # First lets translate the recall list into the positive sample and negative sample list
    pn_sample_list = []
    # Get the first sample
    if recall[0] > 0:
        pn_sample_list.append(1)
    else:
        pn_sample_list.append(0)
    # Get the other samples   
    for i in range(1, len(recall)):
        if recall[i] > recall[i-1]:
            pn_sample_list.append(1)
        else:
            pn_sample_list.append(0)
    # Calculate the number of positive samples and the number of negative samples       
    pn_sample_list = np.array(pn_sample_list)
    positive_samples = np.sum(pn_sample_list)
    negative_samples = len(pn_sample_list)- positive_samples
    # Calculate the TRP and FPR
    TPR = []
    FPR = []
    TP = 0
    FP = 0
    for sample in pn_sample_list:
        if sample == 1:
            TP += 1        
        elif sample == 0:
            FP += 1
            
        TPR.append(TP/positive_samples)
        FPR.append(FP/negative_samples)
    # Calculate the AUC
    AUC = 0
    for i in range(1, len(TPR)):
        AUC += TPR[i]*(FPR[i]- FPR[i-1])
    
    return FPR, TPR, AUC
#########################################################################
def Add_ROC_AUC_to_results(test_coverage_RBFN_results):

    for key, value in test_coverage_RBFN_results.items():
        FPR_list = []
        TPR_list = []
        AUC_list = []
        for recall in value['recalls_list']:
            FPR, TPR, AUC = Recall_to_ROC_AUC(recall)
            FPR_list.append(copy.deepcopy(FPR))
            TPR_list.append(copy.deepcopy(TPR))
            AUC_list.append(AUC)
        # Add those to the results
        value['FPR_list'] = FPR_list
        value['TPR_list'] = TPR_list
        value['AUC_list'] = AUC_list
#########################################################################3
#if __name__ == '__main__':                
#                    
#    test_coverage_RBFN_numerical_results = Calculate_recall_precision(binary = False)  
#
#    reverse_test_coverage_RBFN_numerical_results = Reverse_test(test_coverage_RBFN_numerical_results)  
#    
#    Add_ROC_AUC_to_results(test_coverage_RBFN_numerical_results)
#    Add_ROC_AUC_to_results(reverse_test_coverage_RBFN_numerical_results)
#    # Save the results
#    os.chdir('/home/leo/Documents/Database/Pipeline_New/Codes/Results')
#    with open('test_coverage_RBFN_numerical_results' , 'w') as f:
#        json.dump(test_coverage_RBFN_numerical_results, f, cls = NumpyEncoder)
#    with open('reverse_test_coverage_RBFN_numerical_results', 'w') as f:
#        json.dump(reverse_test_coverage_RBFN_numerical_results, f, cls=NumpyEncoder)
##########################################################################################
'''***********Compare the numerical and binary*************'''
def Draw_ROC():
    os.chdir('/home/leo/Documents/Database/Pipeline_New/Codes/Results')
    with open('test_coverage_RBFN_numerical_results', 'r') as f:
        test_coverage_RBFN_numerical_results = json.load(f)
    with open('test_coverage_RBFN_results', 'r') as f:
        test_coverage_RBFN_results = json.load(f)
    for i in range(1, 4):
        for j in range(1,4):            
            match_type = '('+str(i)+','+str(j)+')'
            save_name ='Numerical'+str(i)+'_'+str(j)
            key = str(i)+'_'+str(j)+'_0_0_1_2_1perchain'
            FPR = test_coverage_RBFN_results[key]['FPR_list'][0]
            TPR = test_coverage_RBFN_results[key]['TPR_list'][0]
            auc = test_coverage_RBFN_results[key]['AUC_list'][0]
            FPR_reverse = test_coverage_RBFN_numerical_results[key]['FPR_list'][0]
            TPR_reverse = test_coverage_RBFN_numerical_results[key]['TPR_list'][0]
            auc_reverse = test_coverage_RBFN_numerical_results[key]['AUC_list'][0]
            plt.figure(figsize = (6, 6))
            plt.plot([0,1], [0,1])
            plt.ylabel('True Positive Rate')
            plt.ylim([0, 1])
            plt.xlabel('False Positive Rate')
            plt.xlim([0,1])
            plt.plot(FPR, TPR, label='Binary:   '+'%.2f'% auc)
            plt.plot(FPR_reverse, TPR_reverse, label='Numerical:   '+'%.2f'% auc_reverse)
            plt.legend(loc=4)
            plt.title('match-type'+match_type+', sample_size =  '+ str(len(FPR)))
#            plt.savefig(save_name+'.png')
            plt.show()

#Draw_ROC()
#Add_ROC_AUC_to_results(test_coverage_RBFN_numerical_results)
#Add_ROC_AUC_to_results(test_coverage_RBFN_results)
#
#for key, value in test_coverage_RBFN_numerical_results.items():
#    print(key, np.average(np.array(value['AUC_list'])))
#
#for key, value in test_coverage_RBFN_results.items():
#    print(key, np.average(np.array(value['AUC_list'])))


