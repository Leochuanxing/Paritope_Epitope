import os
import json
import matplotlib.pyplot as plt
import numpy as np
import copy
import pandas as pd
import random
from matplotlib.ticker import StrMethodFormatter
##########################################################
os.chdir('/home/leo/Documents/Database/Pipeline_New/Cores')

##################################################################3
def Distribution_over_CDR():
    os.chdir('/home/leo/Documents/Database/Pipeline_New/Cores')
    n_for_match_type = []
    paratope_length_distribution = {}
    for key in ['h1', 'h2', 'h3', 'l1', 'l2', 'l3']:
        paratope_length_distribution[key] = []
    for i in range(1, 7):
        l1 = 0; l2=0; l3=0
        h1=0; h2=0; h3=0
        for j in range(1, 7):
            match_type = str(i)+'_'+str(j)
            name = 'training_'+match_type+'_0_0_1_2_1perchain'
            with open(name, 'r') as f:
                cores = json.load(f)
            n_for_match_type.append([match_type, len(cores)])
            for match in cores:
                if match[3][:2] == 'h1':
                    h1 += 1
                elif match[3][:2] == 'h2':
                    h2 += 1
                elif match[3][:2] == 'h3':
                    h3 += 1   
                elif match[3][:2] == 'l1':
                    l1 += 1
                elif match[3][:2] == 'l2':
                    l2 += 1
                elif match[3][:2] == 'l3':
                    l3 += 1
        
        paratope_length_distribution['h1'].append(h1)
        paratope_length_distribution['h2'].append(h2)
        paratope_length_distribution['h3'].append(h3)
        paratope_length_distribution['l1'].append(l1)
        paratope_length_distribution['l2'].append(l2)
        paratope_length_distribution['l3'].append(l3)
        
    core_distribution = {}
    core_distribution['over_match_type'] = n_for_match_type
    core_distribution['over_CDR'] = paratope_length_distribution
    
    return core_distribution
                    
#core_distribution = Distribution_over_CDR()                
#core_distribution                
#os.chdir('/home/leo/Documents/Database/Pipeline_New/Codes/Results')            
#with open('core_distribution', 'w') as f:
#    json.dump(core_distribution, f)
##################################################################
def Distribution_of_aa():
    TripleSingle =  [['TYR', 'Y'], ['LYS', 'K'],['ASP', 'D'], ['ASN', 'N'], ['TRP', 'W'], ['PHE', 'F'], ['GLN', 'Q'],
   ['GLU', 'E'], ['PRO', 'P'], ['GLY', 'G'], ['THR', 'T'],['SER', 'S'], ['ARG', 'R'], ['HIS', 'H'],
   ['LEU', 'L'], ['ILE', 'I'], ['CYS', 'C'], ['ALA', 'A'], ['MET', 'M'], ['VAL', 'V']]
    
    aa_20 = []
    for ts in TripleSingle:
        aa_20.append(ts[0])
        
    distribution_of_aa = {}
    os.chdir('/home/leo/Documents/Database/Pipeline_New/Cores')
    for i in range(1, 4):
        for j in range(1, 4):
            match_type = str(i)+'_'+str(j)
            name = 'training_'+match_type+'_0_0_1_2_1perchain'
            with open(name, 'r') as f:
                cores = json.load(f)
            Ab_aa = []; Ag_aa = []
            for match in cores:
                Ab_aa.extend(match[0])
                Ag_aa.extend(match[1])
                
    # Count the number of each amino acids
            distribution_of_aa[match_type] = {}
            distribution_of_aa[match_type]['Ab_aa'] = []
            distribution_of_aa[match_type]['Ag_aa'] = []
    
            counter_Ab = {}; counter_Ag = {}
            for aa in aa_20:
                counter_Ab[aa]=0; counter_Ag[aa]=0
    
            for a in Ab_aa:
                counter_Ab[a] += 1
            for b in Ag_aa:
                counter_Ag[b] += 1
    
            for key, value in counter_Ab.items():
                distribution_of_aa[match_type]['Ab_aa'].append([key, value])
            for key, value in counter_Ag.items():
                distribution_of_aa[match_type]['Ag_aa'].append([key, value])
            
        
    return distribution_of_aa

#distribution_of_aa = Distribution_of_aa()
#distribution_of_aa

###########################################################
os.chdir('/home/leo/Documents/Database/Pipeline_New/Codes/Results')
#with open('aa_distribution', 'w') as f:
#    json.dump(distribution_of_aa, f) 
           
with open('aa_distribution', 'r') as f:
    aa_distribution = json.load(f)
'''
Draw side by side bar graph
'''
def Draw_aa_distribution(aa_distribution):
    for key, value in aa_distribution.items():
        Ag_aa = value['Ag_aa']
        Ab_aa = value['Ab_aa']
    
        # Change the the frequency into relative frequency
        t_ag = 0; t_ab = 0
        for a in Ag_aa:
            t_ag += a[1]
        for b in Ab_aa:
            t_ab += b[1]
        
        for a in Ag_aa:
            a.append(a[1]/t_ag)
        for b in Ab_aa:
            b.append(b[1]/t_ab)
        
        # Sort Ab_aa
        Ab_aa.sort(key = lambda x:x[1], reverse = True)
        aa_20 = []
        for aa in Ab_aa:
            aa_20.append(aa[0])
        
        Ab_relative = []
        Ag_relative = []
        for aa in aa_20:
            for ab_aa in Ab_aa:
                if ab_aa[0] == aa:
                    Ab_relative.append(ab_aa[1]/t_ab)
            for ag_aa in Ag_aa:
                if ag_aa[0] == aa:
                    Ag_relative.append(ag_aa[1]/t_ag)

        os.chdir('/home/leo/Documents/Database/Pipeline_New/Codes/Results/Data_description')
        x = np.arange(20)
        plt.figure(figsize = (60, 24))
        plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.2f}'))
        plt.bar(x-0.15, Ab_relative, width = 0.3, label = 'Antibody Core Amino Acids')
        plt.bar(x+0.15, Ag_relative, width = 0.3, label = 'Antigen Core Amino Acids')
        plt.xticks(ticks = np.arange(20),labels=aa_20, fontsize=40)
        plt.ylabel('Relative frequency', fontsize = 50)
        plt.rcParams['ytick.labelsize']=40
        plt.legend(loc = 1, prop={'size': 50})
        plt.title(key)
        plt.xlim([-0.5, 20])
        plt.savefig('aa_distribution_'+key+'.eps')
    
#Draw_aa_distribution(aa_distribution)

###################################################################
'''
Get the top match paires of one to one
'''
def Top_one_to_one():
    os.chdir('/home/leo/Documents/Database/Pipeline_New/Cores')
    name = 'training_1_1_0_0_1_2_1perchain'
    with open(name, 'r') as f:
        data = json.load(f)
    # Get the list of the one_to_one
    one = []
    for parepi in data:
        match = []
        match = parepi[0]
        match.extend(parepi[1])
        one.append(copy.deepcopy(match))
    # count the frequencies
    count_one = {}
    for match in one:
        if match[0]+match[1] not in count_one:
            count_one[match[0]+match[1]] = 1
        else:
            count_one[match[0]+match[1]] += 1
    # Generate the one_frequency
    one_frequency = []
    for key, value in count_one.items():
        one_frequency.append([key[:3], key[3:], value])
    # Sort
    one_frequency.sort(key = lambda x:x[2], reverse = True)

    return one_frequency
#one_frequency = Top_one_to_one()  
#len(one_frequency)
#one_frequency  
#os.chdir('/home/leo/Documents/Database/Pipeline_New/Codes/Results') 
#with open('one_to_one_frequency', 'w') as f:
#    json.dump(one_frequency, f)
# Save the one-frequency as csv file
#d = {}
#d['Antibody Amino Acids'] = []
#d['Antigen Amino Acids'] = []
#d['Frequency'] = []
#for mat in one_frequency:
#    d['Antibody Amino Acids'].append(mat[0])
#    d['Antigen Amino Acids'].append(mat[1])
#    d['Frequency'].append(mat[2])
#df = pd.DataFrame(data=d)   
 
####################################################################
def Cores_distribution_over_CDR():
    os.chdir('/home/leo/Documents/Database/Pipeline_New/Codes/Results')            
    with open('core_distribution', 'r') as f:
        core_distribution = json.load(f)
    core_distribution
    
    '''
    Draw Core distribution over CDRs
    '''
    h1 = core_distribution['over_CDR']['h1']
    h2 = core_distribution['over_CDR']['h2']
    h3 = core_distribution['over_CDR']['h3']
    l1 = core_distribution['over_CDR']['l1']
    l2 = core_distribution['over_CDR']['l2']
    l3 = core_distribution['over_CDR']['l3']
    
    #h_matrix = np.array(h1)
    h_matrix = np.stack((h1, h2, h3), axis = 0)
    l_matrix = np.stack((l1, l2, l3), axis = 0)
    h_matrix
    h1
    # Take sum
    h_sum = np.sum(h_matrix, axis = 0)
    l_sum = np.sum(l_matrix, axis = 0)
    # Calculate the relative frequency
    relative_h = h_matrix / h_sum
    relative_l = l_matrix/l_sum
    #relative_h
    # Draw the graph
    x = np.array([0, 1, 2, 3, 4, 5])
    plt.figure(figsize = (40, 24))
    plt.bar(x-0.2, relative_h[2], width = 0.2, label='CDRH3')
    plt.bar(x, relative_h[1], width=0.2, label='CDRH2')
    plt.bar(x+0.2, relative_h[0], width=0.2, label='CDRH1')
    plt.xlim([-0.5, 6])
    plt.legend(loc = 1)
    plt.legend(loc = 1, prop={'size': 40})
    plt.xticks(x, labels = ['1', '2', '3', '4','5', '6'], fontsize = 40)
    plt.rcParams['ytick.labelsize']=40
    plt.xlabel('Length of antibody core amino acids', fontsize = 50)
    plt.ylabel('Relative frequency', fontsize = 50)
    plt.title('Distribution of cores over different CDRHs')
    plt.savefig('distribution_over_CDRH.eps')
    
    x = np.array([0, 1, 2, 3, 4, 5])
    plt.figure(figsize = (40, 24))
    plt.bar(x-0.2, relative_l[2], width = 0.2, label='CDRL3')
    plt.bar(x, relative_l[1], width=0.2, label='CDRL2')
    plt.bar(x+0.2, relative_l[0], width=0.2, label='CDRL1')
    plt.xlim([-0.5, 6])
    plt.legend(loc = 1)
    plt.legend(loc = 1, prop={'size': 40})
    plt.xticks(x, labels = ['1', '2', '3', '4', '5', '6'], fontsize = 40)
    plt.rcParams['ytick.labelsize']=40
    plt.xlabel('Length of antibody core amino acids', fontsize = 40)
    plt.ylabel('Relative frequency', fontsize = 50)
    plt.title('Distribution of cores over different CDRLs')
    plt.savefig('distribution_over_CDRL.eps')
    
#Cores_distribution_over_CDR()    

################################################################
def Core_distribution_over_match_type():
    from mpl_toolkits.mplot3d import axes3d
#    from matplotlib import style
#    style.use('ggplot')
    os.chdir('/home/leo/Documents/Database/Pipeline_New/Codes/Results')
    with open('core_distribution', 'r') as f:
        core_distribution = json.load(f)
    
#    fig = plt.figure(figsize= (40, 24))
    fig = plt.figure(figsize= (10, 6))
    ax1 = fig.add_subplot(111, projection='3d')
    
    
    x3 = []
    y3 = []
    for i in range(1, 7):
        for j in range(1, 7):
            x3.append(i-0.25)
            y3.append(j-0.25)
            
#    x3 = np.array([1,1,1,2,2,2,3,3,3])-0.25
#    y3 = np.array([1,2,3,1,2,3,1,2,3])-0.25
    z3 = np.zeros(len(x3))
    
    dx = np.ones(len(x3))*0.5
    dy = np.ones(len(x3))*0.5
    # Load the heights of the bars
    core_distribution_over_match_type = core_distribution['over_match_type']
    dz = []
    for i in range(1, 7):
        for j in range(1, 7):
            match = str(i)+'_'+str(j)
            for match_type in core_distribution_over_match_type:
                if match_type[0] == match:
                    dz.append(match_type[1])
            
    
    
    ax1.bar3d(x3, y3, z3, dx, dy, dz, color= 'c')
    
    
    ax1.set_xlabel('Length of antibody amino acids')
    ax1.set_ylabel('Length of antigen amino acids')
    ax1.set_zlabel('Number of Cores')
    
    ax1.set_xticks(ticks = [1, 2, 3,4,5,6])
    ax1.set_yticks(ticks = [1,2,3,4,5,6])
    # Save
    os.chdir('/home/leo/Documents/Database/Pipeline_New/Codes/Results')
    plt.savefig('Core_districution_over_Match_type.eps')
    plt.show()
#Core_distribution_over_match_type()
######################################################################
os.chdir('/home/leo/Documents/Database/Pipeline_New/All with peptide 5+ resolution 4A')
with open('training_ac_contact', 'r') as f:
    training_ac_contact = json.load(f)
#
#len(training_ac_contact)
#training_ac_contact
##########################################################################
def The_latest_complex():
    import csv
    os.chdir('/home/leo/Documents/Database/Pipeline_New/Latest')
    pdbid_data = []
    with open('summary.tsv', newline= '') as csvfile:
        total = csv.reader(csvfile, delimiter='\t', quotechar='|')
        for row in total:
            if len(row) >= 10:
                # Make sure the complex is either peptide or protein
                antigen_type = row[5].split('|')
                for at in antigen_type:
                    at.strip()
                if len(row[0]) == 4:
                    if 'protein' in antigen_type or 'peptide' in antigen_type:
                        pdbid_data.append([row[0], row[9], antigen_type])
    # Convert the date to numbers 
    for com in pdbid_data:
        date = com[1]
        date_list = date.split('/')
        if int(date_list[2])>25:
            number_date = int('19'+date_list[2]+date_list[0]+date_list[1])
        else:
            number_date = int('20'+date_list[2]+date_list[0]+date_list[1])
        com[1] = number_date
        
    latest = []
    for date in pdbid_data:
        if date[1] > 20180401 and date[0] not in latest:
            latest.append(date[0])
    
    os.chdir('/home/leo/Documents/Database/Pipeline_New/All with peptide 5+ resolution 4A')
    with open('combined_ids', 'r') as f:
        combined_ids = json.load(f)
    
    real_latest = []
    for pdbid in latest:
        if pdbid not in combined_ids:
            real_latest.append(pdbid)
    return real_latest, pdbid_data
#########################################################################
#import csv
#os.chdir('/home/leo/Documents/Database/Pipeline_New/All with peptide 5+ resolution 4A')
#pdbid_data = []
#with open('summary.tsv', newline= '') as csvfile:
#    total = csv.reader(csvfile, delimiter='\t', quotechar='|')
#    for row in total:
#        if len(row) >= 10:
#            # Make sure the complex is either peptide or protein
#            antigen_type = row[5].split('|')
#            for at in antigen_type:
#                at.strip()
#            if len(row[0]) == 4:
#                if 'protein' in antigen_type or 'peptide' in antigen_type:
#                    pdbid_data.append([row[0], row[9], antigen_type])
## Convert the date to numbers 
#for com in pdbid_data:
#    date = com[1]
#    date_list = date.split('/')
#    if int(date_list[2])>25:
#        number_date = int('19'+date_list[2]+date_list[0]+date_list[1])
#    else:
#        number_date = int('20'+date_list[2]+date_list[0]+date_list[1])
#    com[1] = number_date
#pdbid_data.sort(key = lambda x:x[1])
#pdbid_data[-6:-1]
#with open('testing_ids', 'r') as f:
#    testing_ids = json.load(f)
#testing_id_date = []
#for td in pdbid_data:
#    if td[0] in testing_ids:
#        testing_id_date.append(td)
#testing_id_date
#########################################################################

with open('good_combined_ids', 'r') as f:
    good_combined_ids = json.load(f)
with open('good_matched_ids', 'r') as f:
    good_matched_ids = json.load(f)
#len(good_combined_ids)
#len(good_matched_ids)
#############################################################
#real_latest, pdbid_data = The_latest_complex()
def Generate_combined_matched_ids_the_latest(real_latest, good_combined_ids, good_matched_ids):
    combined_ids_latest = {}
    matched_ids_latest = {}
    for key, value in good_combined_ids.items():
        if key in real_latest:
            combined_ids_latest[key] = value
    for key, value in good_matched_ids.items():
        if key in real_latest:
            matched_ids_latest[key] = value
    return combined_ids_latest, matched_ids_latest

######################################################################################
'''
Draw graphs of the testing and the discrimination testing
'''
def Draw_binary_test_AUC():
    average_AUC={}
    os.chdir('/home/leo/Documents/Database/Pipeline_New/Codes/Results')
    with open('test_results_combined', 'r') as f:
        data = json.load(f)
    
    for i in range(1, 4):
        for j in range(1, 4):
            os.chdir('/home/leo/Documents/Database/Pipeline_New/Codes/Results/Binary_prediction')
            key = str(i)+'_'+str(j)+'_0_0_1_2_1perchain'
            fpr_list = data[key]['FPR_list']
            tpr_list = data[key]['TPR_list']
            average_auc = data[key]['AUC_average']
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
            
            average_AUC[str(i)+'_'+str(j)]= average_auc
    return average_AUC
#average_AUC = Draw_binary_test_AUC()
#average_AUC
#os.chdir('/home/leo/Documents/Database/Pipeline_New/Codes/Results/Binary_prediction')
#with open('average_AUC', 'w') as f:
#    json.dump(average_AUC, f)
         
###########################################################################
def Calculate_FPR_by_TPR(TPR):
    FPR = []
    # The score for the positive is 1 and the score for the negative is 0
    sample = []
    if TPR[0]==0:
        sample.append(0)
    else:
        sample.append(1)
    
    for i in range(1, len(TPR)):
        if TPR[i]-TPR[i-1] == 0:
            sample.append(0)
        else:
            sample.append(1)
        
    # Calculate the total number of positive samples and the total number of negative samples
    total_positive = np.sum(np.array(sample))
    total_negative = len(sample) - total_positive
    # Calculate the FPR
    n=0;p=0
    for i in range(len(sample)):
        if sample[i] == 1:
            p += 1
        if sample[i] == 0:
            n += 1
            
        FPR.append(n/total_negative)
    return FPR

##############################################################################
def Bootstrap_by_TPR(TPR):
    sample_score = []
    if TPR[0]==0:
        sample=-1; score = 0
    else:
        sample = 1; score = 0
    sample_score.append([sample, score])
    
    for i in range(1,len(TPR)):
        if TPR[i]-TPR[i-1] > 0 :
            sample = 1; score = -i
        else:
            sample = -1; score = -i
        sample_score.append([sample, score])
    # Let boot strap
    AUC_list = []
    for r in range(10_000):
        boot_samples = random.choices(sample_score, k=len(sample_score))
        auc = Calculate_AUC(boot_samples)
        AUC_list.append(auc)
    # Calculate the CI_95
    AUC_list.sort()
    CI_95 = [AUC_list[250], AUC_list[9250]]
    return CI_95
    
############################################################################
def Calculate_AUC(sample):
    sample.sort(key = lambda x:x[1], reverse = True)
    # Calculate n_positive and n_negative
    n_positive = 0; n_negative = 0
    for s in sample:
        if s[0] == 1:
            n_positive += 1
        elif s[0] == -1:
            n_negative += 1
    # calculate the tpr and fpr
    neg=0;pos=0; tpr = []; fpr = []
    for s in sample:
        if s[0] == 1:
            pos += 1
        elif s[0] == -1:
            neg += 1
        tpr.append(pos/n_positive)
        fpr.append(neg/n_negative)
    # Calculate the auc
    auc = 0
    for i in range(1, len(tpr)):
        auc += tpr[i]*(fpr[i]-fpr[i-1])        
    return auc
#####################################################################
def Draw_discrimination_AUC():    
    reverse_auc_dict = {}
#    for mode in ['binary', 'latest_binary']:
#    name = 'reverse_test_'+mode
    os.chdir('/home/leo/Documents/Database/Pipeline_New/Codes/Results')
    with open('reverse_test_results', 'r') as f:
        data = json.load(f)
    for i in range(1, 4):
        for j in range(1, 4):
            os.chdir('/home/leo/Documents/Database/Pipeline_New/Codes/Results/Binary_prediction')
            key = str(i)+'_'+str(j)+'_0_0_1_2_1perchain'
            tpr = data[key]['TPR']
            fpr = data[key]['FPR']
            auc = data[key]['AUC']
            CI95 = Bootstrap_by_TPR(tpr)
            # Plot
            plt.figure(figsize = (10, 10))
            plt.plot([0,1], [0,1], linewidth = 2)
            plt.ylabel('True Positive Rate',fontsize = 20)
            plt.ylim([0, 1])
            plt.xlabel('False Positive Rate',fontsize = 20)
            plt.xlim([0,1])
            plt.plot(fpr, tpr, linewidth = 2)
            tic = [0.2, 0.4, 0.6, 0.8, 1]
            lab = ['0.2', '0.4', '0.6', '0.8', '1']
            plt.xticks(tic, labels = lab, fontsize = 15)
            plt.yticks(tic, labels = lab, fontsize = 15)
            plt.rcParams['ytick.labelsize']=10
            plt.rcParams['xtick.labelsize']=10
            plt.title('Match-type('+str(i)+','+str(j)+')', fontsize = 20)
            plt.text(0.4,0.15, 'AUC='+str(round(auc,2))\
                     +'  CI95=['+str(round(CI95[0],2))+','+str(round(CI95[1],2))+']', fontsize = 20)
            plt.savefig('reverse_'+str(i)+'_'+str(j)+'.eps')
            # Calculate the CI95
            print('bootstraping '+ key)
            
            reverse_auc_dict[str(i)+'_'+str(j)]=[auc, CI95]
    return reverse_auc_dict
###########################################################################

reverse_auc_dict_all_negative = Draw_discrimination_AUC() 
#reverse_auc_dict_all_negative
#os.chdir('/home/leo/Documents/Database/Pipeline_New/Codes/Results/Binary_prediction')   
#with open('reverse_auc_dict', 'w') as f:
#    json.dump(reverse_auc_dict_all_negative, f)

############################################################################
'''
****************************************************************************
Analyze the data from the affinity validation
'''


def Affinity_pred():   
    CI95 = {}
    os.chdir('/home/leo/Documents/Database/Pipeline_New/Codes/Results')
    with open('affinity_pred_results', 'r') as f:
        one_to_one_affinity = json.load(f)
    os.chdir('/home/leo/Documents/Database/Pipeline_New/Codes/Results/Affinity_prediction')    
    keys = ['0_1000', '0_0.5', '0.5_1000', '1_1000']
    lengends = [r'$|\Delta\Delta G|>0.0 $', r'$|\Delta\Delta G|<0.5 $',\
               r'$|\Delta\Delta G|>0.5 $', r'$|\Delta\Delta G|>1.0 $']
    plt.figure(figsize = (10, 10))
    plt.plot([0,1], [0,1])
    plt.plot([0,1], [0,1], linewidth = 2)
    plt.ylabel('True Positive Rate', fontsize = 20)
    plt.ylim([0, 1])
    plt.xlabel('False Positive Rate', fontsize = 20)
    plt.xlim([0,1])
    CI = []
    for key, le in zip(keys, lengends):
        fpr = one_to_one_affinity[key]['my_pred_dict']['fpr']
        tpr = one_to_one_affinity[key]['my_pred_dict']['tpr']
        auc = one_to_one_affinity[key]['my_pred_dict']['auc']
        ci = one_to_one_affinity[key]['my_pred_dict']['CI_95']
        CI.append(ci)
        plt.plot(fpr, tpr, label = le+'   '+str(round(auc,2)), linewidth =2)
    tic = [0.2, 0.4, 0.6, 0.8, 1]
    lab = ['0.2', '0.4', '0.6', '0.8', '1']
    plt.xticks(tic, labels = lab, fontsize = 15)
    plt.yticks(tic, labels = lab, fontsize = 15)
    plt.rcParams['ytick.labelsize']=15
    plt.rcParams['xtick.labelsize']=15
    plt.legend(loc = 4, prop={'size': 20})
    plt.title('CM', fontsize = 15)
    plt.savefig('Affinity_pred_'+'CM.eps') 
    CI95['Affinity_pred_'+'CM'] = copy.deepcopy(CI)
    
    methods = ['bASA', 'dDFIRE', 'DFIRE', 'Design_studio', 'FoldX', 'Rosetta', 'Statium']
    for method in methods:
        plt.figure(figsize = (10, 10))
        plt.plot([0,1], [0,1], linewidth = 2)
        plt.ylabel('True Positive Rate', fontsize = 20)
        plt.ylim([0, 1])
        plt.xlabel('False Positive Rate', fontsize = 20)
        plt.xlim([0,1])
        CI = []
        for key, le in zip(keys, lengends):
            fpr = one_to_one_affinity[key]['other_pred_dict'][method]['other_fpr']
            tpr = one_to_one_affinity[key]['other_pred_dict'][method]['other_tpr']
            auc = one_to_one_affinity[key]['other_pred_dict'][method]['other_auc']
            ci = one_to_one_affinity[key]['other_pred_dict'][method]['CI_95']
            CI.append(ci)
            plt.plot(fpr, tpr, label = le+'   '+str(round(auc,2)), linewidth = 2)
        tic = [0.2, 0.4, 0.6, 0.8, 1]
        lab = ['0.2', '0.4', '0.6', '0.8', '1']
        plt.xticks(tic, labels = lab, fontsize = 15)
        plt.yticks(tic, labels = lab, fontsize = 15)
        plt.rcParams['ytick.labelsize']=15
        plt.rcParams['xtick.labelsize']=15
        plt.legend(loc = 4, prop={'size': 20})
        plt.title(method, fontsize = 15)
        plt.savefig('Affinity_pred_'+method+'.eps')
        CI95['Affinity_pred_'+method] = copy.deepcopy(CI)
    return CI95
#CI95 = Affinity_pred()
#os.chdir('/home/leo/Documents/Database/Pipeline_New/Codes/Results/Affinity_prediction') 
#with open('CI95', 'w') as f:
#    json.dump(CI95, f)
##################################################################################
def Affinity_pred_on_my_data():
    os.chdir('/home/leo/Documents/Database/Pipeline_New/Codes/Results')
    with open('prediction_on_my_data', 'r') as f: 
        prediction_on_my_data = json.load(f)
    keys = ['0_1000', '0_0.5', '0.5_1000', '1_1000']
    legends = [r'$|\Delta\Delta G|>0.0 $', r'$|\Delta\Delta G|<0.5 $',\
           r'$|\Delta\Delta G|>0.5 $', r'$|\Delta\Delta G|>1.0 $']
    os.chdir('/home/leo/Documents/Database/Pipeline_New/Codes/Results/Affinity_prediction')
    plt.figure(figsize = (30, 30))
    plt.plot([0,1], [0,1], linewidth = 2)
    plt.ylabel('True Positive Rate', fontsize = 60)
    plt.ylim([0, 1])
    plt.xlabel('False Positive Rate', fontsize = 60)
    plt.xlim([0,1])
    CI95_on_my_data = []
    for key, le in zip(keys, legends):
        fpr = prediction_on_my_data[key]['fpr']
        tpr = prediction_on_my_data[key]['tpr']
        auc = prediction_on_my_data[key]['auc']
        ci =  prediction_on_my_data[key]['CI_95']
        CI95_on_my_data.append([auc,ci])
        plt.plot(fpr, tpr, label = le+'   '+str(round(auc,2)), linewidth = 3)
    tic = [0.2, 0.4, 0.6, 0.8, 1]
    lab = ['0.2', '0.4', '0.6', '0.8', '1']
    plt.xticks(tic, labels = lab, fontsize = 40)
    plt.yticks(tic, labels = lab, fontsize = 40)
    plt.rcParams['ytick.labelsize']=40
    plt.rcParams['xtick.labelsize']=40
    plt.legend(loc = 4, prop={'size': 50})
    plt.title('MC on my data', fontsize = 30)
    plt.savefig('MC_on_my_data.eps')
        
    return CI95_on_my_data
        
#CI95_on_my_data = Affinity_pred_on_my_data()    
#with open('CI95_on_my_data', 'w') as f:
#    json.dump(CI95_on_my_data, f)
################################################################################
'''
Some basic numbers
'''
















