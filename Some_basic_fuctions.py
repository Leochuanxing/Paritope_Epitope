'''THIS FILE IS ABOUT SOME FREQUENTLY USED FUNCTIONS'''
import random
#######################################################################################
'''
AUC_TPR_FPR:
    This function is to calculate AUC, TPR, FPR
Input:
    pred: a list, gives the predicted values
    observed; a list, gives the observed values
    cut: a float, if the observed values <= cut, it is defined as a negative sample.
                  If the observed values > cut, they are defined as the positive samples
                  
    REMARK: 
        1, The negative samples are defined as the samples with the observed values 
    
                   *****less than or equal to cut*****
                   
        2. The pred and the observed should corresponding with each other in order.
        
Output:
    AUC: area under the curve
    TPR: true positive rate, a list
    FPR: false positive rate, a list
'''

def AUC_TPR_FPR(pred, observed_values, cut):
     # Calculate the number of positive observations
    positive_total = 0 ; negative_total = 0
    for observed in observed_values:
        if observed > cut :
            positive_total += 1
        elif observed <= cut:
            negative_total += 1
    if negative_total == 0:
        print('None negative samples in Calculating AUC')
        return None, None, None
    if positive_total == 0:
        print('None positive samples in Calculating AUC')
        return None, None, None
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
        if match[1] > cut:
            n_positive += 1
        elif match[1] <= cut:
            n_negative += 1
            
        TPR.append(n_positive/positive_total)
        FPR.append(n_negative/negative_total)
    
    # Calculate the AUC
    AUC = 0
    for i in range(1, len(TPR)):
        AUC += TPR[i] * (FPR[i] - FPR[i-1])
    
    return AUC, TPR, FPR
######################################################################################
'''
Input: 
    pred, observed, cut are as described in the AUC_TPR_FPR function.
    iteration: a positive integer, gives the number of bootstrap
Output:
    AUC_list:
        gives a list of AUC values, of which the length equals iteration.
'''
def Bootstrap_AUC(pred, observed, cut, iteration):
    n = len(pred)
    # Sample with replacement
    AUC_list = []
    for i in range(iteration):
        bootstrap_ind = random.choices(list(range(n)), k=n)
        bootstrap_pred = [pred[i] for i in bootstrap_ind]
        bootstrap_observed = [observed[i] for i in bootstrap_ind]
        # Calculate AUC
        AUC, _, _ = AUC_TPR_FPR(bootstrap_pred, bootstrap_observed, cut)
        # Load the data
        AUC_list.append(AUC)
        
    return AUC_list






































