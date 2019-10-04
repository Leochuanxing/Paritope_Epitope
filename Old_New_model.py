'''OLD MODEL AND NEW MODEL'''
import os
import json
os.chdir('/home/leo/Documents/Database/Data_Code_Publish/Codes')
from Predict_affinity import Predict_affinity, Analyze_resutls
#############################################################################
os.chdir('/home/leo/Documents/Database/Pipeline_New/Complexes/Results')
with open('one_to_one_affinity', 'r') as f:
    one_to_one_affinity = json.load(f)
one_to_one_affinity.keys()
one_to_one_affinity['affinity_results']['results_dict'].keys()    
one_to_one_affinity.keys()
len(one_to_one_affinity['1_1000']['others']['fpr'])
one_to_one_affinity['0.5_1000']['auc_all']
# Extract the information
results = one_to_one_affinity['affinity_results']['results_dict']['single_True_mode_single_moving_True_binary_True']['0']
workable_old = []
mut_sets=[]; DDG = []; mut_id = []
for res in results:
    one_mut_set = []; one_workable = []
    for i in range(1, len(res)-1):
        if res[i] != []:
            one_mut_set.append([res[i][0][0], res[i][1][0]]) 
    if one_mut_set != []:
        for mut in one_mut_set:
            one_workable.append([mut[0],mut[1], 1, res[0][2][2], res[0][3]]) 
        mut_sets.append(one_mut_set)
        DDG.append(res[0][2][2])
        mut_id.append(res[0][3])
        workable_old.append(one_workable)
#workable_old[:10]
''' [...
[[wt_Ab_Ag_pair], [mut_Ab_Ag_pair], contact_number, DDG, mut_id],
...]'''
#working_d = '/home/leo/Documents/Database/Data_Code_Publish/Codes/Results'
#pred_res_old = Predict_affinity(workable_old, working_d, binary = True)
#selected_cut_DDG, AUC, TPR, FPR, correct_ratio = \
#            Analyze_resutls(pred_res_old, cut_DDG_lower=1, cut_DDG_upper=100)
#
#AUC_old = AUC
#correct_ratio_old = correct_ratio
#AUC_old 
#correct_ratio_old
#AUC
#len(TPR)
##pred_res_old
#########################################################################################
#'''CHECH THE EXTRACTED MUTATIONS'''
#saving_d = '/home/leo/Documents/Database/Data_Code_Publish/Codes/Results'
#os.chdir(saving_d) 
#with open('workable_output', 'r') as f:
#    workable_n = json.load(f)       
#
#workable_new = workable_n['multiple_WithinRange_True']  
#len(workable)
#len(workable_new)            
#ids = []
#for i in workable:
#    ids.append(i[0][-1])  
#
#len(ids)  
#for ind, m in enumerate(ids):
#    if m == 210:
#        print(ind)
#i = 31
#for mu in workable_new:
#    if mu != [] and mu[0][-1] == ids[i]:
#        print(mu)
#for mu_old in workable:
#    if mu_old != [] and mu_old[0][-1] == ids[i]:
#        print(mu_old)
#            
#results[:10]
###############################################################################
#os.chdir('/home/leo/Documents/Database/Data_Code_Publish/Structures')
#with open('combined_ids', 'r') as f:
#    combined_ids = json.load(f)            
#            
#combined_ids['1bj1']           

            
            
            
            
            
            
            
            
            
            