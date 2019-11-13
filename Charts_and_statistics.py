import numpy as np
import copy
import pandas as pd
import random
import math
import json
import os
from scipy.stats import norm
from scipy.stats import ks_2samp
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
    
results_d = '/home/leo/Documents/Database/Data_Code_Publish/Codes/Results'
save_d = '/home/leo/Documents/Database/Data_Code_Publish/Codes/Results/Graphs'
#Draw_all(results_d, save_d)
    
'''***********************************************************************'''    
''' SECTION, SECTION SECTION SECTION SECTION SECTION SECTION SECTION SECTION '''
'''***********************************************************************'''
    
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
                Ag_tuple.append((aa, freq, round(freq/Ag_n, 3), Ag_n))
            Ag_tuple.sort(key=lambda x:x[1], reverse=True)
            core_aa_freq[str(i)+'_'+str(j)]['Ag_aa'] = Ag_tuple[:]    

    # Store the above results in a dictionary
    AA_freq = {}
    AA_freq['Ag_aa'] = Ag_aa_tuple_list
    AA_freq['Ab_aa'] = Ab_aa_tuple_list
    AA_freq['CDR_aa_freq'] = CDR_aa_freq
    AA_freq['core_aa_freq'] = core_aa_freq             
                    
    return AA_freq


def Two_sample_z_proportion(n1, p1, n2, p2, test = 'p1>p2'):
    # Calculat the pooled proportion
    p_pool = (n1*p1 + n2*p2)/(n1+n2)
    variance = p_pool*(1-p_pool)*(1/n1 + 1/n2)
    # Calculate the test statistic z
    z = (p1-p2)/(variance**0.5)
    # Calculate p_value
    if test == 'p1>p2':
        p_value = 1 - norm.cdf(z)
    elif test == 'p1<p2':
        p_value = norm.cdf(z)
    elif test == 'p1 != p2':
        p_value = 2*min(norm.cdf(z), 1-norm.cdf(z))
    else:
        print('test is not proper!')
        return
        
    return p_value

def Sub_Top_5_aa_distribution(aa_match_type, base_line, top = 5):
    top_n = []
    for k in range(top):
        aa_info_match_type = list(aa_match_type[k])
        for aa_base_line in base_line:
            if aa_base_line[0] == aa_info_match_type[0]:
                aa_info_base_line = aa_base_line
                break
        # Calculate the p values for Ab-aa
        n1 = aa_info_match_type[3]
        p1 = aa_info_match_type[2]
        n2 = aa_info_base_line[3]
        p2 = aa_info_base_line[2]
        p_value = Two_sample_z_proportion(n1,p1,n2,p2,test= 'p1>p2')
        
        aa_info_match_type.append(p_value)
        top_n.append(aa_info_match_type)
        
    return top_n

def Top_5_aa_distribution(AA_freq):
    # Compare the aa distribution in cores with aa distribution in CDRs
    Ag_aa = AA_freq['Ag_aa']
    CDR_aa_freq = AA_freq['CDR_aa_freq']
    core_aa_freq = AA_freq['core_aa_freq']
    # Add up the aa in CDR_aa_freq
    sum_aa_CDR = {}; total = 0
    for CDR, aa in CDR_aa_freq.items():
        for aa_freq in aa:
            if aa_freq[0] not in sum_aa_CDR:
                sum_aa_CDR[aa_freq[0]] = aa_freq[1]
                total += aa_freq[1]
            else:
                sum_aa_CDR[aa_freq[0]] += aa_freq[1]
                total += aa_freq[1]
    # Calculate the relative frequency of aa in CDRs
    aa_CDR = []
    for aa, freq in sum_aa_CDR.items():
        aa_CDR.append((aa, freq, round(freq/total, 5), total))
    aa_CDR.sort(key=lambda x:x[1])
    # Calculate the p values of the top 5 aa of different match type 
    top_5_aa = {}; all_aa = {}
    for i in range(1, 4):
        for j in range(1, 4):
            match_type = str(i)+'_'+str(j)
            top_5_aa[match_type] = {}
            all_aa[match_type] = {}
            
            Ab_aa_match_type = core_aa_freq[match_type]['Ab_aa']
            Ag_aa_match_type = core_aa_freq[match_type]['Ag_aa']
            
            all_aa[match_type]['Ab'] = Sub_Top_5_aa_distribution(Ab_aa_match_type, aa_CDR, top = 20)
            all_aa[match_type]['Ag'] = Sub_Top_5_aa_distribution(Ag_aa_match_type, Ag_aa, top = 20)
            
            top_5_aa[match_type]['Ab'] = all_aa[match_type]['Ab'][:5]
            top_5_aa[match_type]['Ag'] = all_aa[match_type]['Ag'][:5]  
            
            all_aa[match_type]['Ab'].sort(key=lambda x:x[4])
            all_aa[match_type]['Ag'].sort(key=lambda x:x[4])
    
    return top_5_aa, all_aa

'''
top_5_aa:
    For each match type, top_5_aa gives the statistics of comparing the top 5 most
    frequent amino acids with the baseline distribution. The baseline distribution of 
    Antibody amino acids is the distribution of the amino acides from the CDRs. The 
    baseline distribution of the antigen amino acids is the distribution of the amino acids
    from the whole antigen.
all_aa:
    It gives all the statistics of all the 20 amino acids. The baselines are the same as
    in the top 5 aa.
AA_freq:
    Gives the frequencies and relative frequencies of the amino acids in cores, CDRs and 
    the whole proterin.
'''

def AA_distribution_statistics():
    AA_freq = AA_relative_freq()
    top_5_aa, all_aa = Top_5_aa_distribution(AA_freq)
    # Load all the above results to one dictionary
    aa_distribution_stats ={}
    aa_distribution_stats['AA_freq'] = AA_freq 
    aa_distribution_stats['top_5_aa'] = top_5_aa
    aa_distribution_stats['all_aa'] = all_aa
    # Save the results
    os.chdir('/home/leo/Documents/Database/Data_Code_Publish/Codes/Results')
    with open('aa_distribution_statistics', 'w') as f:
        json.dump(aa_distribution_stats, f)

        


'''***********************************************************************'''    
''' SECTION, SECTION SECTION SECTION SECTION SECTION SECTION SECTION SECTION '''
'''***********************************************************************'''


def Count_cores_over_CDRs():
    l_range = [[23, 40], [49, 63], [89, 110]]
    h_range = [[25, 37], [50, 71], [99, 129]]

    l_total = 0; h_total = 0
    for l, h in zip(l_range, h_range):
        l_total += l[1] - l[0] +1
        h_total += h[1] - h[0] +1
        
    l_prop = []; h_prop = []    
    for l, h in zip(l_range, h_range):
        l_prop.append((l[1]-l[0]+1)/l_total)
        h_prop.append((h[1]-h[0]+1)/h_total)
        
    
    l_exp_ob = {}; h_exp_ob = {}
    os.chdir('/home/leo/Documents/Database/Data_Code_Publish/Cores/Positive_cores')
    core_count_CDRs = {}
    for i in range(1,4):
        for j in range(1, 4):
            match_type = str(i)+'_'+str(j)
    
            test_name = 'testing_'+match_type+'_0_0_1_2_1perchain'
            train_name = 'training_'+match_type+'_0_0_1_2_1perchain'
            
            with open(test_name,'r') as f:
                testing = json.load(f)
            with open(train_name,'r') as f:
                training= json.load(f)
        
            cores_match_type = testing[:]
            cores_match_type.extend(training)
            
            h1 = 0; h2 = 0; h3 = 0; l1 = 0; l2 = 0; l3 = 0
            for fcdn in cores_match_type:
                if fcdn[3][:2] == 'h1':
                    h1 += 1
                elif fcdn[3][:2] == 'h2':
                    h2 += 1 
                elif fcdn[3][:2] == 'h3':
                    h3 += 1 
                elif fcdn[3][:2] == 'l1':
                    l1 += 1 
                elif fcdn[3][:2] == 'l2':
                    l2 += 1 
                elif fcdn[3][:2] == 'l3':
                    l3 += 1 
                    
            core_count_CDRs[match_type] = [(h1, 'h1'),(h2, 'h2'),(h3, 'h3'),(l1, 'l1'),\
                           (l2, 'l2'),(l3, 'l3')]
            
            l_exp_ob[match_type] = []
            h_exp_ob[match_type] = []
            
            l_total = l1 + l2 + l3
            l_exp_ob[match_type].append([l_total*l_prop[0], l_total*l_prop[1], l_total*l_prop[2]])
            l_exp_ob[match_type].append([l1, l2, l3])
            
            h_total = h1 + h2 + h3
            h_exp_ob[match_type].append([h_total*h_prop[0], h_total*h_prop[1], h_total*h_prop[2]])
            h_exp_ob[match_type].append([h1, h2, h3])

            
    return core_count_CDRs, l_exp_ob, h_exp_ob, l_prop, h_prop


def One_p_z_test(l_exp_ob, l_prop):
    exp_prop = l_prop[2]
    pvs = {}
    for i in range(1,4):
        for j in range(1, 4):
            match_type = str(i)+'_'+str(j)
            ob = l_exp_ob[match_type][1]
            cdr_total = ob[0]+ob[1]+ob[2]
            ob_prop = ob[2]/cdr_total
            z = ob_prop-exp_prop
            z /= (exp_prop*(1-exp_prop)/cdr_total)**0.5
            p_value = 1 - norm.cdf(z)
            
            pvs[match_type] = p_value
    return pvs


def Top_5_percent_cores():
    os.chdir('/home/leo/Documents/Database/Data_Code_Publish/Cores/Positive_cores')
    top_5_percent_cores = {}
    for i in range(1, 4):
        for j in range(1,4):
            match_type = str(i)+'_'+str(j)
            test_name = 'testing_'+match_type+'_0_0_1_2_1perchain'
            train_name = 'training_'+match_type+'_0_0_1_2_1perchain'
            
            
            with open(test_name, 'r') as f:
                test = json.load(f)
            with open(train_name, 'r') as f:
                train = json.load(f)
            
            # Extract the top 5 percent cores, frequency and relative frequency
            cores_list = test[:]
            cores_list.extend(train)
            
            cores = []
            for fcdn in cores_list:
                AA=''
                for Ab_aa in fcdn[0]:
                    AA += Ab_aa
                for Ag_aa in fcdn[1]:
                    AA += Ag_aa
                cores.append(AA)
                
            cores_set = set(cores)
            
            cores_count = []
            for c in cores_set:
                cores_count.append([c, cores.count(c)])
            cores_count.sort(key=lambda x:x[1], reverse=True)
            
            # Append the relative frequency
            total = len(cores_list)
            for cc in cores_count:
                cc.append(round(cc[1]/total, 5))
                
            temp = cores_count[:math.ceil(0.05*len(cores_set))]
            # Change the cores into normal form
            container = []
            for core_freq in temp:
                core_verb = core_freq[0]
#                l = len(core_verb)
                Ab_verb = core_verb[:i*3]
                Ag_verb = core_verb[i*3:]
                Ab_aa = []; Ag_aa = []
                for k in range(i):
                    Ab_aa.append(Ab_verb[3*k:3*(k+1)])
                for k in range(j):
                    Ag_aa.append(Ag_verb[3*k:3*(k+1)])
                container.append([Ab_aa, Ag_aa, core_freq[1], core_freq[2]])
            # Load to the dictionary
            top_5_percent_cores[match_type]=container[:]
                
    return top_5_percent_cores


def Core_over_CDR_statistics():
    core_count_CDRs, l_exp_ob, h_exp_ob, l_prop, h_prop = Count_cores_over_CDRs()
    chains = ['heavy_chain', 'light_chain']
    p_values_CDR3_more_cores = {}
    for chain in chains:
        if chain == 'heavy_chain':
            p_values_CDR3_more_cores[chain]=One_p_z_test(h_exp_ob, h_prop)
        elif chain == 'light_chain':
            p_values_CDR3_more_cores[chain] = One_p_z_test(l_exp_ob, l_prop)
    # Load all the results into a dicitonary
    core_over_CDR_statistics = {}
    core_over_CDR_statistics['p_values_CDR3_more_cores'] = p_values_CDR3_more_cores
    core_over_CDR_statistics['core_count_CDRs'] = core_count_CDRs    

    top_5_percent_cores = Top_5_percent_cores()  
    core_over_CDR_statistics['top_5_percent_cores'] = top_5_percent_cores    

    # Save the results
    os.chdir('/home/leo/Documents/Database/Data_Code_Publish/Codes/Results')
    with open('core_over_CDR_statistics', 'w') as f:
        json.dump(core_over_CDR_statistics, f)
    
    
def Statistics():
    AA_distribution_statistics()
    Core_over_CDR_statistics()
    
#Statistics()


'''***********************************************************************'''    
''' SECTION, SECTION SECTION SECTION SECTION SECTION SECTION SECTION SECTION '''
'''***********************************************************************'''

def Core_distribution_over_match_type():
#    from mpl_toolkits.mplot3d import axes3d  
#    fig = plt.figure(figsize= (40, 24))
    fig = plt.figure(figsize= (10, 6))
    ax1 = fig.add_subplot(111, projection='3d')
    
    
    x3 = []
    y3 = []
    for i in range(1, 5):
        for j in range(1, 5):
            x3.append(i-0.25)
            y3.append(j-0.25)

    z3 = np.zeros(len(x3))
    
    dx = np.ones(len(x3))*0.5
    dy = np.ones(len(x3))*0.5
    # Load the heights of the bars
    os.chdir('/home/leo/Documents/Database/Data_Code_Publish/Cores/Positive_cores')
    dz = []; nCores_match_type = []
    for i in range(1, 5):
        for j in range(1, 5):
            match = str(i)+'_'+str(j)
            test_name = 'testing_'+match+'_0_0_1_2_1perchain'
            train_name = 'training_'+match+'_0_0_1_2_1perchain'
            with open(test_name, 'r') as f:
                test = json.load(f)
            with open(train_name, 'r') as f:
                train = json.load(f)
            dz.append(len(test)+len(train))
            nCores_match_type.append((match, len(test)+len(train)))
    
    ax1.bar3d(x3, y3, z3, dx, dy, dz, color= 'c')
    
    
    ax1.set_xlabel('Length of antibody amino acids')
    ax1.set_ylabel('Length of antigen amino acids')
    ax1.set_zlabel('Number of Cores')
    
    ax1.set_xticks(ticks = [1, 2, 3,4,5])
    ax1.set_yticks(ticks = [1,2,3,4,5])
    # Save
    os.chdir('/home/leo/Documents/Database/Data_Code_Publish/Codes/Results/Graphs')
    plt.savefig('Core_districution_over_Match_type.eps')
    plt.show()
    
    os.chdir('/home/leo/Documents/Database/Data_Code_Publish/Codes/Results')
    with open('nCores_match_type', 'w') as f:
        json.dump(nCores_match_type, f)

#Core_distribution_over_match_type()

'''***********************************************************************'''    
''' SECTION, SECTION SECTION SECTION SECTION SECTION SECTION SECTION SECTION '''
'''***********************************************************************'''

def Compare_Ab_aa_Ag_aa_in_cores():

    os.chdir('/home/leo/Documents/Database/Data_Code_Publish/Codes/Results')
    with open('aa_distribution_statistics', 'r') as f:
        aa_distribution_stats = json.load(f)
    Ab_aa = {}; Ag_aa = {}
    for i in range(1, 4):
        for j in range(1, 4):
            key = str(i)+'_'+str(j)
            for aa in aa_distribution_stats['AA_freq ']['core_aa_freq'][key]['Ab_aa'][:20]:
                if aa[0] not in Ab_aa:
                    Ab_aa[aa[0]] = aa[1]
                else:
                    Ab_aa[aa[0]] += aa[1]
            for aa in aa_distribution_stats['AA_freq ']['core_aa_freq'][key]['Ag_aa'][:20]:
                if aa[0] not in Ag_aa:
                    Ag_aa[aa[0]] = aa[1]
                else:
                    Ag_aa[aa[0]] += aa[1]
    # Calculate total                
    Ab_total = 0; Ag_total = 0
    for aa, n in Ab_aa.items():
        Ab_total += n
    for aa, n in Ag_aa.items():
        Ag_total += n
    # Calculate relative frequency
    Ab_Ag_relfreq = []
    for aa, n in Ab_aa.items():
        Ab_Ag_relfreq.append([aa, round(n/Ab_total, 3), round(Ag_aa[aa]/Ag_total, 3), n, Ag_aa[aa]])
        
    # Sort accroding to Ab
    Ab_Ag_relfreq.sort(key = lambda x:x[1], reverse = True)
    # Load the results to the aa distribution dicitonary
    aa_distribution_stats['Ab_Ag_in_cores_felfreq'] = Ab_Ag_relfreq
    # Save the results 
    
    os.chdir('/home/leo/Documents/Database/Data_Code_Publish/Codes/Results')
    with open('aa_distribution_statistics', 'w') as f:
        json.dump(aa_distribution_stats, f)
    
    # Draw graph and save
        
    Ab_relative = []; Ag_relative = []; aa_20 = []
    for Ab_Ag in Ab_Ag_relfreq:
        aa_20.append(Ab_Ag[0])
        Ab_relative.append(Ab_Ag[1])
        Ag_relative.append(Ab_Ag[2])

        
    x = np.arange(20)
    plt.figure(figsize = (15, 6))
    #plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.2f}'))
    plt.bar(x-0.15, Ab_relative, width = 0.3, label = 'Antibody Core Amino Acids')
    plt.bar(x+0.15, Ag_relative, width = 0.3, label = 'Antigen Core Amino Acids')
    plt.xticks(ticks = np.arange(20),labels=aa_20, fontsize=15)
    plt.ylabel('Relative frequency', fontsize = 30)
    plt.rcParams['ytick.labelsize']=20
    plt.legend(loc = 1, prop={'size': 25})
    #plt.title('Antibody amino acids and Antibody amino acids in cores', y=1.05, fontsize = 40)
    plt.xlim([-0.5, 20])
    os.chdir('/home/leo/Documents/Database/Data_Code_Publish/Codes/Results/Graphs')
    plt.savefig('aa_distribution_overall.eps')
    plt.show()

#Compare_Ab_aa_Ag_aa_in_cores()
































