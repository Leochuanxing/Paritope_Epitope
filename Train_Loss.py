'''
############THIS FILE IS THE LOSS FUNCTION AND THE TRAINING METHODS################
This file may be constantly used in machine learning. The function will be imported is

                ****************Train*******************   
                
The input of this funciton is train_para, it should contain the following:
   
train_para['observed']:
    It gives the observed values. In dimension (n, k)
train_para['design_matrix']: 
    The design matrix, after the kernel transformations. In dimension (n, m)
train_para['coefficients']:
    The inner product between the above design matrix and this coefficents will 
    be the final results of the model. In dimension (m ,k)
train_para['reg']: 
    The regularization coefficient, to control the overfitting. This could be a hyperparameter
train_para['loss_type']:
    It takes the value from {'Sigmoid', 'Softmax', 'SVM', 'SumSquares'}, deciding the loss 
    function will be used in the training process
train_para['method']:
    It takes the value from {'BFGS', 'GD'}, where BFGS is the Quasi Newton method with Ramijo
    back track line search, and GD means gradient descent.
train_para['step_size'];
    If the method takes the value of GD, step_size must be assigned. It gives the length of each in 
    Gradient descent.
train_para['n_iterations']:
    If the method takes the value of GD, n_iterations must be assigned. It gives 
    the number of iterations before the GD process ends.
 
    
            **************************************************
            ****************************************************
            THE FOLLOW VALUES CAN BE COPIED AND PASTED FOR CONVENIENCE
            
train_para['observed']
train_para['design_matrix']
train_para['coefficients']
train_para['reg']
train_para['loss_type']
train_para['method']
train_para['step_size']
train_para['n_iterations']
    
'''

import numpy as np

'''#################################################################################
****************************************************************************************
                CALCULATE THE LOSS AND THE GRADIENTS

Suppose: design_matrix is of dimension (m,n)                

Loss_Sigmoid: Use the sigmoid as the loss function
Inputs:
    design_matrix: as above
    labels: an array corresponding to the rows of the design_matrix,
            with values of either 0 or 1. labels.shape = (m, 1)
    coefficients: an array contains the coeffients, coefficients.shape = (n, 1)
    reg: a float, the coefficient of regularization
Outputs:
    loss: float
    grad_coefficients: an array containing all the gradients corresopnding to the coefficients
'''
def Loss_Sigmoid(design_matrix, labels, coefficients, reg):
    
    nrow, ncol = design_matrix.shape
    
    logit = design_matrix.dot(coefficients)
    prob = 1/(1+np.exp(-logit))
    loss = np.average(- np.log(prob) * labels - (1 - labels) * np.log(1 - prob))
    # plus the regularization
    loss += reg * np.sum(coefficients * coefficients)
    
    # Calculate the gradient from the first part of loss
    grad_logit = prob - labels
    grad_coefficients = (design_matrix.T).dot(grad_logit)
    grad_coefficients /= nrow
    # Calculate the gradient from the regularizatio part
    grad_coefficients += 2 * reg * coefficients
    
    # return the above results
    return loss, grad_coefficients
'''
Loss_Softmax: this function applies to the case when the output classes are more than 2
Input:
    design_matrix: as above
    labels: a matrix of dimension (m, k), where k is the number of classes. The label of
            each sample is a vector of dimenion (1, k), and the values are either 0 or 1, with
            1 indicate the correct category.
    coefficients: a matrix of dimension (n*k, 1). For the convenience of latter usage,  we don't use the shape(n, k).
                    When (n,k) is reshaped into (n*k,1), we stack column by column.
    reg: as above
Output:
    similar as above
    
THE FLAW OF SOFTMAX: IF THE DIFFERENCES OF THE LOGIT VALUES ARE TWO BIG, THE LOSS FUNCTION MAY BE TOO BIG!
'''
def Loss_Softmax(design_matrix, labels, coefficients, reg):
    
    nrow, ncol = design_matrix.shape
    # Reshape the coefficients
    coefficients = coefficients.reshape((-1, ncol)).T
    
    Wx = design_matrix.dot(coefficients)
    # Make sure the elements in Wx is not too big or too small
    Wx -= np.max(Wx, axis = 1, keepdims = True)
    # Calculate the probabilities
    exp = np.exp(Wx)
    prob = exp / np.sum(exp, axis = 1, keepdims = True)
    
    log_prob = np.log(prob)

    # Calculate  the loss
    loss = np.sum(- log_prob * labels)/nrow
    loss += reg * np.sum(coefficients * coefficients)
    
    # Calculate the gradients
    grad_Wx = prob - labels
    grad_coefficients = (design_matrix.T).dot(grad_Wx)
    grad_coefficients /= nrow
    
    grad_coefficients += 2 * reg * coefficients
    
    grad_coefficients = grad_coefficients.T.reshape((-1, 1))
    
    return loss, grad_coefficients

'''
'''

def Loss_SVM(design_matrix, observed, coefficients, reg):
     
    nrow, ncol = design_matrix.shape
    # Reshape the coefficients
    coefficients = coefficients.reshape((-1, ncol)).T
    # Calculate the loss
    ii = np.zeros((observed.shape[1], observed.shape[1])) + 1
    Wx = design_matrix.dot(coefficients)
    s1 = Wx + 1
    obs = observed * Wx
    obsii = obs.dot(ii)
    ad = s1 - obsii
    d = ad * (1-observed)
    ind = (d>0)
    sd = d * ind
    loss = np.sum(sd)
    loss += reg * np.sum(coefficients * coefficients)
    
    # Calculate the gradients
    grad_d = ind
    grad_ad = grad_d * (1-observed)
    grad_s1 = grad_ad
    grad_obsii = - grad_ad
    grad_Wx = grad_s1
    grad_obs = grad_obsii.dot(ii)
    grad_Wx += observed * grad_obs
    grad_coeff = (design_matrix.T).dot(grad_Wx)
    
    grad_coeff += 2 * reg * coefficients
    # Reshape the gradient
    grad_coeff = grad_coeff.T.reshape((-1, 1))
    return loss, grad_coeff


'''
Loss_SumSquares: the loss is measured by the sum of the sequares of the difference 
                between the predicted values and the observed values
Inputs:
    observed: an array of shape (m,1). with each element a float.
    design_matrix, coefficients and reg are the same as above
Outputs:
    the same as above
'''    
def Loss_SumSquares(design_matrix, observed, coefficients, reg):
    
    nrow, ncol = design_matrix.shape
    
    # Calculate the loss
    pred = design_matrix.dot(coefficients)
    loss = np.average((pred - observed) * (pred - observed))
    
    loss += reg * np.sum(coefficients * coefficients)
    
    # Calculate the gradient
    
    grad_coefficients = (design_matrix.T).dot(2 * (pred - observed))
    grad_coefficients /= nrow
    grad_coefficients += 2 * reg * coefficients
    
    return loss, grad_coefficients
'''
Integrate the above functions into one function for convenient usage.
Input:
    train_para: a dictionary, contains all the needed values to train a model.
'''
def Loss(train_para):
    design_matrix = train_para['design_matrix'] 
    observed = train_para['observed']
    reg = train_para['reg']
    coefficients = train_para['coefficients']
    loss_type = train_para['loss_type']
    if loss_type == 'Sigmoid':
        loss, grad_coefficients = Loss_Sigmoid(design_matrix, observed, coefficients, reg)
    elif loss_type == 'Softmax':
        loss, grad_coefficients = Loss_Softmax(design_matrix, observed, coefficients, reg)
    elif loss_type == 'SumSquares':
        loss, grad_coefficients = Loss_SumSquares(design_matrix, observed, coefficients, reg)
    elif loss_type == 'SVM':
        loss, grad_coefficients = Loss_SVM(design_matrix, observed, coefficients, reg)
    return loss, grad_coefficients


'''#################################################################################'''
'''
Train_GD: train the model by using gradient descent
Inputs:
    gd_train_para: a dictionary, contains
        reg: coefficients of regularization
        setp_size: a float
        loss_type: string, gives types of different loss functions
        design_matrix:
        observed: observed values, the format of which decides the type of loss functions
Outputs:  
    coefficients: a matrix of shape (design_matrix.shape[0], observed.shape[1])
                    the values of the coefficients after n_iterations training      
'''
def Train_GD(train_para):
    # Take out the parameters
#    design_matrix = train_para['design_matrix'] 
#    observed = train_para['observed']
#    reg = train_para['reg']
    step_size = train_para['step_size']
#    loss_type = gd_train_para['loss_type']
    n_iterations = train_para['n_iterations']    
    
    for i in range(n_iterations):
        loss, grad_coefficients = Loss(train_para)
        # update the coefficients
        train_para['coefficients'] -= step_size * grad_coefficients
        '''Do we print the loss'''
        if i % 100 == 0:
            print(round(loss, 6))
            
    return train_para, loss

'''
The parameters in Train_RBFN_BFGS are explained in related textbook about quasi Newton
method and Ramijo research.
'''
def Train_RBFN_BFGS(train_para, rho=0.8, c = 1e-4, termination = 1e-2):   
    
    nrow, _ = np.shape(train_para['coefficients']) 
    max_design = np.max(np.abs(train_para['design_matrix']))       
    # Create an iteration counter
    n_iteration = 0
    # BFGS algorithm
    loss, grad_coeff = Loss(train_para)
    # Initiate H. This H should not be large in case it may destroy the Loss function
    H = np.eye(nrow)
    H *= 1/(np.max(np.abs(grad_coeff)) * max_design)
    ternination_square = termination**2
    grad_square = ternination_square + 1
    while grad_square >= ternination_square:          
        # keep a record of this grad_square for monitoring the efficiency of this process
        n_iteration += 1
      
        p = - H.dot(grad_coeff)        
        # There should be both old and new coefficients in the train_para
        train_para['coefficients_old'] = train_para['coefficients']
        train_para['coefficients'] = p + train_para['coefficients_old']
        # Calculate the loss and gradient
        new_loss, new_grad_coeff = Loss(train_para)        
        # Ramijo Back-tracking
        while new_loss > loss + c * (grad_coeff.T).dot(p):
            p *= rho
            train_para['coefficients'] = p + train_para['coefficients_old']            
            new_loss, new_grad_coeff = Loss(train_para)        
        # update H
        s = p
        y = new_grad_coeff - grad_coeff
        r = (y.T).dot(s)
        I = np.eye(nrow)
        if r != 0:
            r = 1/r            
            H = (I - r*s.dot(y.T)).dot(H).dot(I - r*y.dot(s.T)) + r*s.dot(s.T)# Can be accelerated
        else:
            H = np.diag(np.random.uniform(0.5, 1, nrow))# try to eliminate the periodic dead loop
            H *= 1/(np.max(np.abs(new_grad_coeff))*max_design)# Make sure H is not too large
        # Update loss, grad_square and paramter
        loss = new_loss
        grad_coeff = new_grad_coeff
        grad_square = new_grad_coeff.T.dot(new_grad_coeff)            
        # print some values to monitor the training process 
        if n_iteration % 100 == 0:
            print('loss  ', loss, '    ','grad_square   ', grad_square)
            n_iteration = 0        
        
    return train_para, loss

'''
Combine the above two training functions
'''
def Train(train_para, rho=0.8, c = 1e-4, termination = 1e-2):
    method = train_para['method']
    if method == 'GD':
        train_para, loss = Train_GD(train_para)
    elif method == 'BFGS':
        train_para, loss = Train_RBFN_BFGS(train_para, rho=0.8, c = 1e-4, termination = 1e-2)
    return train_para, loss
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    