import numpy as np
import copy
import pandas as pd
import random
import math
import json
import os
#from matplotlib.ticker import StrMethodFormatter
import matplotlib.pyplot as plt
#######################################################################
# Draw auroc curves 
'''
Draw graphs of the testing and the discrimination testing
'''
def My_AUROC_figure(match_type, fpr, tpr, average_auc):
    i = match_type[0]
    j = match_type[1]
    
    roc = plt.figure(figsize = (10, 10))
    plt.plot([0,1], [0,1])
    plt.ylabel('True Positive Rate', fontsize = 20)
    plt.ylim([0, 1])
    plt.xlabel('False Positive Rate', fontsize = 20)
    plt.xlim([0,1])
    tic = [0.2, 0.4, 0.6, 0.8, 1]
    lab = ['0.2', '0.4', '0.6', '0.8', '1']
    plt.xticks(tic, labels = lab, fontsize = 15)
    plt.yticks(tic, labels = lab, fontsize = 15)
    plt.rcParams['ytick.labelsize']=10
    plt.rcParams['xtick.labelsize']=10
#    for fpr, tpr in zip(fpr_list, tpr_list):
    plt.plot(fpr, tpr)
    plt.title('Match-type('+str(i)+','+str(j)+')', fontsize = 20)
    plt.text(0.6, 0.2, 'Average AUC: ' + str(round(average_auc, 2)), fontsize = 20)
    
    return roc

    
def Draw_binary_test_AUC(results_d, save_d):
#    average_AUC={}
    os.chdir(results_d)
    with open('test_discrimination', 'r') as f:
        test_discrimination = json.load(f)

    binary = ['binary', 'numerical']
    for model in binary:        
        data = test_discrimination['cross_'+model+'_Gaussian_']
        for i in range(1, 5):
            for j in range(1, 5):
                os.chdir(save_d)
                key = str(i)+'_'+str(j)+'_0_0_1_2_1perchain'
                fpr_list = data[key]['FPR_list']
                tpr_list = data[key]['TPR_list']
                average_auc = np.average(np.array(data[key]['AUC_list']))
        #        auc_list = data[key]['AUC_list']
                # Draw the plots
                plt.figure(figsize = (10, 10))
                plt.plot([0,1], [0,1])
                plt.ylabel('True Positive Rate', fontsize = 20)
                plt.ylim([0, 1])
                plt.xlabel('False Positive Rate', fontsize = 20)
                plt.xlim([0,1])
                tic = [0.2, 0.4, 0.6, 0.8, 1]
                lab = ['0.2', '0.4', '0.6', '0.8', '1']
                plt.xticks(tic, labels = lab, fontsize = 15)
                plt.yticks(tic, labels = lab, fontsize = 15)
                plt.rcParams['ytick.labelsize']=10
                plt.rcParams['xtick.labelsize']=10
                for fpr, tpr in zip(fpr_list, tpr_list):
                    plt.plot(fpr, tpr)
                plt.title('Match-type('+str(i)+','+str(j)+')', fontsize = 20)
                plt.text(0.6, 0.2, 'Average AUC: ' + str(round(average_auc, 2)), fontsize = 20)
                plt.savefig(str(i)+'_'+str(j)+'.eps')
                

def Draw_integrated_graphs(results_d, save_d):
    os.chdir(results_d)
    with open('integrated_discrimination', 'r') as f:
        data = json.load(f)

    for i in range(1, 4):
        for j in range(1, 4):
            os.chdir(save_d)
            key = 'testing_'+str(i)+'_'+str(j)+'_0_0_1_2_1perchain'
            fpr_list = data[key]['FPR_list']
            tpr_list = data[key]['TPR_list']
            average_auc = np.average(np.array(data[key]['AUC_list']))
    #        auc_list = data[key]['AUC_list']
            # Draw the plots
            plt.figure(figsize = (10, 10))
            plt.plot([0,1], [0,1])
            plt.ylabel('True Positive Rate', fontsize = 20)
            plt.ylim([0, 1])
            plt.xlabel('False Positive Rate', fontsize = 20)
            plt.xlim([0,1])
            tic = [0.2, 0.4, 0.6, 0.8, 1]
            lab = ['0.2', '0.4', '0.6', '0.8', '1']
            plt.xticks(tic, labels = lab, fontsize = 15)
            plt.yticks(tic, labels = lab, fontsize = 15)
            plt.rcParams['ytick.labelsize']=10
            plt.rcParams['xtick.labelsize']=10
            for fpr, tpr in zip(fpr_list, tpr_list):
                plt.plot(fpr, tpr)
            plt.title('Match-type('+str(i)+','+str(j)+')', fontsize = 20)
            plt.text(0.6, 0.2, 'Average AUC: ' + str(round(average_auc, 2)), fontsize = 20)
            plt.savefig('integrated_'+str(i)+'_'+str(j)+'.eps')
            
            
def Draw_test_reverse_AUC(results_d, save_d):
#    average_AUC={}
    os.chdir(results_d)
    with open('test_reverse', 'r') as f:
        test_reverse = json.load(f)

    binary = ['binary', 'numerical']
    for model in binary:        
        data = test_reverse['cross_'+model+'_Gaussian_']
        for i in range(1, 5):
            for j in range(1, 5):
                os.chdir(save_d)
                key = str(i)+'_'+str(j)+'_0_0_1_2_1perchain'
                fpr = data[key]['FPR']
                tpr = data[key]['TPR']
                average_auc = np.average(np.array(data[key]['AUC']))
        #        auc_list = data[key]['AUC_list']

                roc = My_AUROC_figure((i,j), fpr, tpr, average_auc)
                plt.show(roc)
                roc.savefig('reverse_'+str(i)+'_'+str(j)+'.eps')
                
                
def Draw_integrated_reverse_AUC(results_d, save_d):
#    average_AUC={}
    os.chdir(results_d)
    with open('integrated_reverse', 'r') as f:
        data = json.load(f)
        
    for i in range(1, 4):
        for j in range(1, 4):
            os.chdir(save_d)
            key = 'testing_'+str(i)+'_'+str(j)+'_0_0_1_2_1perchain'
            fpr = data[key]['FPR']
            tpr = data[key]['TPR']
            average_auc = np.average(np.array(data[key]['AUC']))

            roc = My_AUROC_figure((i,j), fpr, tpr, average_auc)
            plt.show(roc)
            roc.savefig('integrated_reverse_'+str(i)+'_'+str(j)+'.eps')
            
            
def Draw_affinity(results_d, save_d):
    os.chdir('/home/leo/Documents/Database/Data_Code_Publish/Codes/Results')
    with open('pred_workable_formated', 'r') as f:
        pred_workable_formated = json.load(f)
    
    os.chdir(save_d)
    lengends = [r'$|\Delta\Delta G|<0.5 $', r'$|\Delta\Delta G|>0.0 $', \
               r'$|\Delta\Delta G|>0.5 $', r'$|\Delta\Delta G|>1.0 $']
    
    WithinRange = ['one_WithinRange_True', 'one_WithinRange_False']
    methods = ['CN', 'bASA', 'dDfire', 'dfire', 'discovery_studio', 'rosetta', 'statium', 'foldX']
    for method in methods:
        for Range in WithinRange:
            plt.figure(figsize = (10, 10))
            plt.plot([0,1], [0,1], linewidth = 2)
            plt.ylabel('True Positive Rate', fontsize = 20)
            plt.ylim([0, 1])
            plt.xlabel('False Positive Rate', fontsize = 20)
            plt.xlim([0,1])
            
            fpr_list = pred_workable_formated[Range][method]['FPR_list']
            tpr_list = pred_workable_formated[Range][method]['TPR_list']
            auc_list = pred_workable_formated[Range][method]['AUC_list']
#            ci_list = pred_workable_formated[Range][method]['BootStrap']

#            CI = []
            for i in range(4):
                fpr = fpr_list[i]
                tpr = tpr_list[i]
                auc = auc_list[i]
#                ci = ci_list[i]
                plt.plot(fpr, tpr, label =lengends[i] +'   '+str(round(auc,2)), linewidth = 2)
            tic = [0.2, 0.4, 0.6, 0.8, 1]
            lab = ['0.2', '0.4', '0.6', '0.8', '1']
            plt.xticks(tic, labels = lab, fontsize = 15)
            plt.yticks(tic, labels = lab, fontsize = 15)
            plt.rcParams['ytick.labelsize']=15
            plt.rcParams['xtick.labelsize']=15
            plt.legend(loc = 4, prop={'size': 20})
            plt.title(method, y = 1.03,fontsize = 50)
            plt.savefig(Range+'_'+method+'.eps')

def Draw_all(results_d, save_d):
    Draw_binary_test_AUC(results_d, save_d)
    Draw_integrated_graphs(results_d, save_d)
    Draw_test_reverse_AUC(results_d, save_d)
    Draw_integrated_reverse_AUC(results_d, save_d)
    Draw_affinity(results_d, save_d)
    
'''***********************************************************************'''    
'''***********************************************************************'''
'''***********************************************************************'''
    
#os.chdir('/home/leo/Documents/Database/Data_Code_Publish/Codes/Results')
#with open('pred_workable_formated', 'r') as f:
#    pred_workable_formated = json.load(f)
#pred_workable_formated.keys()
#pred_workable_formated['cut_range']
#pred_workable_formated['iteration']
#pred_workable_formated['one_WithinRange_True']['CN'].keys()
#pred_workable_formated['one_WithinRange_True']['CN']['BootStrap']
#len(pred_workable_formated['one_WithinRange_True']['CN']['TPR_list'])
#pred_workable_formated['one_WithinRange_True'].keys()

results_d = '/home/leo/Documents/Database/Data_Code_Publish/Codes/Results'
save_d = '/home/leo/Documents/Database/Data_Code_Publish/Codes/Results/Graphs'
#Draw_all(results_d, save_d)
'''
This output of AA_relative_freq should be the relative frequency and the frequency 
of different amino acids for antibody, antigen, cdr, and cores
'''
def AA_relative_freq():
    # Calculate the frequency and relative frequency for antibody and antigens
    # Use the complexes after the redundency elimination
    os.chdir('/home/leo/Documents/Database/Data_Code_Publish/Structures')
    with open('testing_ac_contact', 'r') as f:
        testing_ac_contact = json.load(f)
    with open('training_ac_contact', 'r') as f:
        training_ac_contact = json.load(f)
    with open('sequence', 'r') as f:
        sequence = json.load(f)
    with open('combined_ids', 'r') as f:
        combined_ids = json.load(f)
    
    complexes = list(testing_ac_contact.keys())
    complexes.extend(list(training_ac_contact))
    
    # Calculate the frequency of the heavy and light chains
    Ab_aa_freq = {}; pdb_list = []; n = 0
    for pdbHL in complexes:
        pdb = pdbHL[:4]; HL = pdbHL[4:5]
        pdb_list.append(pdb)
        seq = sequence[pdb][HL]
        for aa in seq:
            if aa  in Ab_aa_freq:
                Ab_aa_freq[aa] += 1
                n += 1
            else:
                Ab_aa_freq[aa] = 1
                n += 1
    # Calculate the relative frequency
    Ab_aa_tuple_list = []
    for key, v in Ab_aa_freq.items():
        Ab_aa_tuple_list.append((key, v, round(v/n, 3), n))
    Ab_aa_tuple_list.sort(key=lambda x:x[1], reverse=True)

    # Calculate the aa frequencies of the antigens
    Ag_aa_freq = {}; n = 0
    pdb_list = list(set(pdb_list))
    for pdb in pdb_list:
        antigen_ids = combined_ids[pdb][2]
        for antigen in antigen_ids:
            anti_seq = sequence[pdb][antigen]
            for aa in anti_seq:
                if aa in Ag_aa_freq:
                    Ag_aa_freq[aa] += 1
                    n += 1
                else:
                    Ag_aa_freq[aa] = 1
                    n += 1
        # Calculate the relative frequency
    Ag_aa_tuple_list = []
    for key, v in Ag_aa_freq.items():
        Ag_aa_tuple_list.append((key, v, round(v/n, 3), n))
    Ag_aa_tuple_list.sort(key=lambda x:x[1], reverse=True)
    
    # Calculate the aa frequencies in the CDR
    l_range = [[23, 40], [49, 63], [89, 110]]
    h_range = [[25, 37], [50, 71], [99, 129]]
    
    l_inds = []; h_inds = []
    for l_ran, h_ran in zip(l_range, h_range):
        l_inds.extend(range(l_ran[0], l_ran[1]+1))
        h_inds.extend(range(h_ran[0], h_ran[1] + 1))
        
    CDR_aa_freq = {}
    CDR_aa_freq['CDRH1'] = {}
    CDR_aa_freq['CDRH2'] = {}
    CDR_aa_freq['CDRH3'] = {}
    CDR_aa_freq['CDRL1'] = {}
    CDR_aa_freq['CDRL2'] = {}
    CDR_aa_freq['CDRL3'] = {}
    
    for pdbHL in complexes:
        pdb = pdbHL[:4]; HL = pdbHL[4:5]; hl = pdbHL[5]
        seq = sequence[pdb][HL]
        if hl == 'h':
            for i in range(len(seq)):
                for j in [0, 1, 2]:
                    key = 'CDRH'+str(j+1)
                    if i in range(h_range[j][0], h_range[j][1] + 1):
                        if seq[i] in CDR_aa_freq[key]:
                            CDR_aa_freq[key][seq[i]] += 1
                        else:
                            CDR_aa_freq[key][seq[i]] = 1
        elif hl == 'l':
            for i in range(len(seq)):
                for j in [0, 1, 2]:
                    key = 'CDRL'+str(j+1)
                    if i in range(l_range[j][0], l_range[j][1] + 1):
                        if seq[i] in CDR_aa_freq[key]:
                            CDR_aa_freq[key][seq[i]] += 1
                        else:
                            CDR_aa_freq[key][seq[i]] = 1
    # Calculate the relative frequencies
    for CDR, aa_freq in CDR_aa_freq.items():
        temp_tuple = []; n = 0
        for aa, freq in aa_freq.items():
            n += freq
        for aa, freq in aa_freq.items():
            temp_tuple.append((aa, freq, round(freq/n, 3), n))
        temp_tuple.sort(key = lambda x:x[1], reverse=True)   
        CDR_aa_freq[CDR] = temp_tuple[:]
                            
    # Calculate the aa frequencies in cores of different match-types
    core_aa_freq = {}
    os.chdir('/home/leo/Documents/Database/Data_Code_Publish/Cores/Positive_cores')
    for i in range(1, 4):
        for j in range(1, 4):
            core_aa_freq[str(i)+'_'+str(j)] = {}
            core_aa_freq[str(i)+'_'+str(j)]['Ab_aa'] = {}
            core_aa_freq[str(i)+'_'+str(j)]['Ag_aa'] = {}
            
            test_key = 'testing_'+str(i)+'_' + str(j)+'_0_0_1_2_1perchain'
            train_key = 'training_'+str(i)+'_' + str(j)+'_0_0_1_2_1perchain'
            with open(test_key, 'r') as f:
                test_cores = json.load(f)
            with open(train_key, 'r') as f:
                train_cores = json.load(f)
                
            cores = []
            cores.extend(test_cores)
            cores.extend(train_cores)
            
            Ab_n = 0; Ag_n = 0
            for match in cores:
                for Ab_aa in match[0]:
                    if Ab_aa in core_aa_freq[str(i)+'_'+str(j)]['Ab_aa']:
                        core_aa_freq[str(i)+'_'+str(j)]['Ab_aa'][Ab_aa] += 1
                        Ab_n += 1
                    else:
                        core_aa_freq[str(i)+'_'+str(j)]['Ab_aa'][Ab_aa] = 1
                        Ab_n += 1
                        
                for Ag_aa in match[1]:
                    if Ag_aa in core_aa_freq[str(i)+'_'+str(j)]['Ag_aa']:
                        core_aa_freq[str(i)+'_'+str(j)]['Ag_aa'][Ag_aa] += 1
                        Ag_n += 1
                    else:
                        core_aa_freq[str(i)+'_'+str(j)]['Ag_aa'][Ag_aa] = 1  
                        Ag_n += 1
            # Calculate the relative frequencies
            Ag_tuple = []; Ab_tuple = []
            for aa, freq in core_aa_freq[str(i)+'_'+str(j)]['Ab_aa'].items():
                Ab_tuple.append((aa, freq, round(freq/Ab_n, 3), Ab_n))
            Ab_tuple.sort(key = lambda x:x[1], reverse=True)
            core_aa_freq[str(i)+'_'+str(j)]['Ab_aa'] = Ab_tuple[:]
            
            for aa, freq in core_aa_freq[str(i)+'_'+str(j)]['Ag_aa'].items():
                Ag_tuple.append((aa, freq, round(freq/Ab_n, 3), Ag_n))
            Ag_tuple.sort(key=lambda x:x[1], reverse=True)
            core_aa_freq[str(i)+'_'+str(j)]['Ag_aa'] = Ag_tuple[:]    

    # Store the above results in a dictionary
    AA_freq = {}
    AA_freq['Ag_aa'] = Ag_aa_tuple_list
    AA_freq['Ab_aa'] = Ab_aa_tuple_list
    AA_freq['CDR_aa_freq'] = CDR_aa_freq
    AA_freq['core_aa_freq'] = core_aa_freq             
                    
    return AA_freq

'''*****************************************************************************'''


