pimport os
import json
import matplotlib.pyplot as plt
import numpy as np
import copy
import pandas as pd
import random
import math
from matplotlib.ticker import StrMethodFormatter
##########################################################
#os.chdir('/home/leo/Documents/Database/Pipeline_New/Cores')
'''Combine the testing sets and stor them in a new folder'''
#for i in range(1, 4):
#    for j in range(1, 4):
#        match_type = str(i)+'_'+str(j)
#        for n in range(10):
##        training_name = 'training_'+match_type+'_0_0_1_2_1perchain'
#            testing_name = 'testing_'+match_type+'_0_0_1_2_1perchain_negative'
#            directory1 = '/home/leo/Documents/Database/Pipeline_New/Negative_Cores/Sample_'+str(n)
#            os.chdir(directory1)
#    #        with open(training_name, 'r') as f:
#    #            training_set = json.load(f)
#            with open(testing_name, 'r') as f:
#                testing_set = json.load(f)
#            
#            directory2 = '/home/leo/Documents/Database/Pipeline_New/Latest/Negative_cores/sample_'+str(n)
#            os.chdir(directory2)
#            with open(testing_name, 'r') as f:
#                testing_set_latest = json.load(f)
#            testing_set.extend(testing_set_latest)
#            # Save the training and the testing set in the new folder
#            os.chdir('/home/leo/Documents/Database/Pipeline_New/Cores/NegativeCores')
#    #        with open('training_'+match_type, 'w') as f:
#    #            json.dump(training_set, f)
#            with open('testing_'+match_type+'_'+str(n), 'w') as f:
#                json.dump(testing_set, f)
##            print('training  ', len(training_set))
#        print('testing1  ', len(testing_set_latest))
#        print('testing2  ', len(testing_set))
#os.chdir('/home/leo/Documents/Database/Pipeline_New/Codes/Results') 
#with open('test_results', 'r') as f:
#    results = json.load(f)     
#publication_res = {}
#for key, value in results.items():
#    new_key = key[:3]
#    publication_res[new_key] = {}
#    publication_res[new_key]['coeff'] = value['coeff']
#    publication_res[new_key]['centers_selected'] = value['centers_selected']
#publication_res.keys()
#publication_res['2_1'].keys()
#len(publication_res['2_1']['coeff'])
#len(publication_res['2_1']['centers_selected'])
#publication_res['2_1']['centers_selected']
# save the results in binary
#os.chdir('/home/leo/Documents/Database/Pipeline_New/Codes/Results/Binary_prediction')
#with open('Centers and coefficients', 'w') as f:
#    json.dump(publication_res, f)
##################################################################3
def Distribution_over_CDR(core_d):
    os.chdir(core_d)
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

##################################################################
def Distribution_of_aa(core_d):
    TripleSingle =  [['TYR', 'Y'], ['LYS', 'K'],['ASP', 'D'], ['ASN', 'N'], ['TRP', 'W'], ['PHE', 'F'], ['GLN', 'Q'],
   ['GLU', 'E'], ['PRO', 'P'], ['GLY', 'G'], ['THR', 'T'],['SER', 'S'], ['ARG', 'R'], ['HIS', 'H'],
   ['LEU', 'L'], ['ILE', 'I'], ['CYS', 'C'], ['ALA', 'A'], ['MET', 'M'], ['VAL', 'V']]
    
    aa_20 = []
    for ts in TripleSingle:
        aa_20.append(ts[0])
        
    distribution_of_aa = {}
    os.chdir(core_d)
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


###########################################################

           

'''
Draw side by side bar graph
'''
def Draw_aa_distribution(results_d):
    os.chdir(results_d)
    with open('aa_distribution', 'r') as f:
        aa_distribution = json.load(f)
        
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


        x = np.arange(20)
        plt.figure(figsize = (15, 6))
        plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.2f}'))
        plt.bar(x-0.15, Ab_relative, width = 0.3, label = 'Antibody Core Amino Acids')
        plt.bar(x+0.15, Ag_relative, width = 0.3, label = 'Antigen Core Amino Acids')
        plt.xticks(ticks = np.arange(20),labels=aa_20, fontsize=15)
        plt.ylabel('Relative frequency', fontsize = 30)
        plt.rcParams['ytick.labelsize']=20
        plt.legend(loc = 1, prop={'size': 25})
        plt.title('Match-type('+key[0]+','+key[-1]+')', y=1.05, fontsize = 40)
        plt.xlim([-0.5, 20])
        plt.savefig('aa_distribution_'+key+'.eps')
    

###################################################################
'''
Get the top match paires of one to one
'''
def Top_one_to_one(core_d, results_d):
    os.chdir(core_d)
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
    # Save the results
    with open('one_to_one_frequency', 'w') as f:
        json.dump(one_frequency, f)
    # Save to csv
    os.chdir(results_d)
    d = {}
    d['Antibody Amino Acids'] = []
    d['Antigen Amino Acids'] = []
    d['Frequency'] = []
    for mat in one_frequency:
        d['Antibody Amino Acids'].append(mat[0])
        d['Antigen Amino Acids'].append(mat[1])
        d['Frequency'].append(mat[2])
    df = pd.DataFrame(data=d)   
    df.to_csv('one_one_frequency.csv')





####################################################################
def Cores_distribution_over_CDR(results_d):
    os.chdir(results_d)            
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
def Core_distribution_over_match_type(results_d):
    from mpl_toolkits.mplot3d import axes3d

    os.chdir(results_d)
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
    os.chdir('/home/leo/Documents/Database/Pipeline_New/Complexes/Results')
    plt.savefig('Core_districution_over_Match_type.eps')
    plt.show()

######################################################################

#def The_latest_complex():
#    import csv
#    os.chdir('/home/leo/Documents/Database/Pipeline_New/Latest')
#    pdbid_data = []
#    with open('summary.tsv', newline= '') as csvfile:
#        total = csv.reader(csvfile, delimiter='\t', quotechar='|')
#        for row in total:
#            if len(row) >= 10:
#                # Make sure the complex is either peptide or protein
#                antigen_type = row[5].split('|')
#                for at in antigen_type:
#                    at.strip()
#                if len(row[0]) == 4:
#                    if 'protein' in antigen_type or 'peptide' in antigen_type:
#                        pdbid_data.append([row[0], row[9], antigen_type])
#    # Convert the date to numbers 
#    for com in pdbid_data:
#        date = com[1]
#        date_list = date.split('/')
#        if int(date_list[2])>25:
#            number_date = int('19'+date_list[2]+date_list[0]+date_list[1])
#        else:
#            number_date = int('20'+date_list[2]+date_list[0]+date_list[1])
#        com[1] = number_date
#        
#    latest = []
#    for date in pdbid_data:
#        if date[1] > 20180401 and date[0] not in latest:
#            latest.append(date[0])
#    
#    os.chdir('/home/leo/Documents/Database/Pipeline_New/All with peptide 5+ resolution 4A')
#    with open('combined_ids', 'r') as f:
#        combined_ids = json.load(f)
#    
#    real_latest = []
#    for pdbid in latest:
#        if pdbid not in combined_ids:
#            real_latest.append(pdbid)
#    return real_latest, pdbid_data
#########################################################################

#########################################################################
#
#with open('good_combined_ids', 'r') as f:
#    good_combined_ids = json.load(f)
#with open('good_matched_ids', 'r') as f:
#    good_matched_ids = json.load(f)
##len(good_combined_ids)
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
def Draw_binary_test_AUC(results_d):
    average_AUC={}
    os.chdir(results_d)
    with open('test_results', 'r') as f:
        data = json.load(f)
    
    for i in range(1, 4):
        for j in range(1, 4):
            os.chdir(results_d)
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

    with open('average_AUC', 'w') as f:
        json.dump(average_AUC, f)

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
def Draw_discrimination_AUC(results_d):    
    reverse_auc_dict = {}
#    for mode in ['binary', 'latest_binary']:
#    name = 'reverse_test_'+mode
    os.chdir(results_d)
    with open('reverse_test_results', 'r') as f:
        data = json.load(f)
    for i in range(1, 4):
        for j in range(1, 4):
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
            
    with open('reverse_auc_dict', 'w') as f:
        json.dump(reverse_auc_dict, f)

###########################################################################

#reverse_auc_dict = Draw_discrimination_AUC() 
#reverse_auc_dict_all_negative
#os.chdir('/home/leo/Documents/Database/Pipeline_New/Complexes/Results')   


############################################################################
'''
****************************************************************************
Analyze the data from the affinity validation
'''
def Fomalize_results(results_d):
    os.chdir(results_d)
    with open('one_to_one_affinity', 'r') as f:
        one_to_one_affinity = json.load(f)
    one_to_one_affinity.keys()
    one_to_one_affinity['0_1000'].keys()
    one_to_one_affinity['0_1000']['my_pred_dict'].keys()
    one_to_one_affinity['0_1000']['other_pred_dict'].keys()
    one_to_one_affinity['0_1000']['other_pred_dict']['bASA'].keys()
    type(one_to_one_affinity['affinity_results'])
    one_to_one_affinity['affinity_results'].keys()
    type(one_to_one_affinity['affinity_results']['results_dict'])
    one_to_one_affinity['affinity_results']['results_dict'].keys()
    # Format the one_to_one_affinity
    one_to_one ={}
    keys = ['0_1000', '0_0.5', '0.5_1000', '1_1000']
    methods = ['bASA', 'dDFIRE', 'DFIRE', 'Design_studio', 'FoldX', 'Rosetta', 'Statium']
    for key in keys:
        one_to_one[key]={}
        one_to_one[key]['CM']={}
        for k, v in one_to_one_affinity[key]['my_pred_dict'].items():
            one_to_one[key]['CM'][k] = v
    for key in keys:
        for md in methods:
            one_to_one[key][md] = {}
            for ke, va in one_to_one_affinity[key]['other_pred_dict'][md].items():
                if len(ke)==9:
                    one_to_one[key][md][ke[-3:]] = va
                elif len(ke) == 14:
                    one_to_one[key][md][ke[-8:]] = va
                else:
                    one_to_one[key][md][ke] = va
    one_to_one.keys()
    one_to_one['0.5_1000'].keys()        
    one_to_one['0.5_1000']['DFIRE'].keys()  
    one_to_one['0.5_1000']['bASA']['auc']
    one_to_one['0_1000']['dDFIRE']['fpr']
    one_to_one['0_0.5']['Design_studio']['tpr']
    one_to_one['0_1000']['FoldX']['CI_95']
    one_to_one['0_1000']['Rosetta']['selected']
    one_to_one['affinity_results'] = one_to_one_affinity['affinity_results']
    
    with open('one_to_one', 'w') as f:
        json.dump(one_to_one, f)




def Affinity_pred(results_d):   
    CI95 = {}
    os.chdir(results_d)
    with open('one_to_one', 'r') as f:
        one_to_one = json.load(f)
    
    keys = ['0_1000', '0_0.5', '0.5_1000', '1_1000']
    lengends = [r'$|\Delta\Delta G|>0.0 $', r'$|\Delta\Delta G|<0.5 $',\
               r'$|\Delta\Delta G|>0.5 $', r'$|\Delta\Delta G|>1.0 $']
    
    methods = ['bASA', 'dDFIRE', 'DFIRE', 'Design_studio', 'FoldX', 'Rosetta', 'Statium', 'CM']
    for method in methods:
        plt.figure(figsize = (10, 10))
        plt.plot([0,1], [0,1], linewidth = 2)
        plt.ylabel('True Positive Rate', fontsize = 20)
        plt.ylim([0, 1])
        plt.xlabel('False Positive Rate', fontsize = 20)
        plt.xlim([0,1])
        CI = []
        for key, le in zip(keys, lengends):
            fpr = one_to_one[key][method]['fpr']
            tpr = one_to_one[key][method]['tpr']
            auc = one_to_one[key][method]['auc']
            ci = one_to_one[key][method]['CI_95']
            CI.append(ci)
            plt.plot(fpr, tpr, label = le+'   '+str(round(auc,2)), linewidth = 2)
        tic = [0.2, 0.4, 0.6, 0.8, 1]
        lab = ['0.2', '0.4', '0.6', '0.8', '1']
        plt.xticks(tic, labels = lab, fontsize = 15)
        plt.yticks(tic, labels = lab, fontsize = 15)
        plt.rcParams['ytick.labelsize']=15
        plt.rcParams['xtick.labelsize']=15
        plt.legend(loc = 4, prop={'size': 20})
        plt.title(method, y = 1.03,fontsize = 50)
        plt.savefig('Affinity_pred_'+method+'.eps')
        CI95['Affinity_pred_'+method] = copy.deepcopy(CI)
    with open('CI95', 'w') as f:
        json.dump(CI95, f)

##################################################################################
def Affinity_pred_on_my_data():
    os.chdir('/home/leo/Documents/Database/Pipeline_New/Complexes/Results')
    with open('prediction_on_my_data', 'r') as f: 
        prediction_on_my_data = json.load(f)
    keys = ['0_1000', '0_0.5', '0.5_1000', '1_1000']
    legends = [r'$|\Delta\Delta G|>0.0 $', r'$|\Delta\Delta G|<0.5 $',\
           r'$|\Delta\Delta G|>0.5 $', r'$|\Delta\Delta G|>1.0 $']
    os.chdir('/home/leo/Documents/Database/Pipeline_New/Complexes/Results')
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
'''Analyze the distribution of different amino acids'''
def Top_concentration(core_d, results_d, percentage=0.05):
    from mpl_toolkits.mplot3d import axes3d
    top_concentration = []
    os.chdir(core_d)
    x3 = []; y3 = []; dz = []
    for i in range(1, 4):
        for j in range(1, 4):
            x3.append(i-0.25)
            y3.append(j-0.25)
            
            name = 'training_'+str(i)+'_'+str(j)+'_0_0_1_2_1perchain'
            with open(name, 'r') as f:
                core = json.load(f)
            # Find the top six matches
            matches = []; match_set = []
            for parepi in core:
                matches.append([parepi[0], parepi[1]])
                if [parepi[0], parepi[1]] not in match_set:
                    match_set.append([parepi[0], parepi[1]])
            # count the number
            for match in match_set:
                n = 0
                for parepi in core:
                    if parepi[0]==match[0] and parepi[1]==match[1]:
                        n+=1
                match.append(n)
            match_set.sort(key = lambda x:x[2], reverse = True)
            n = math.floor(len(match_set)*percentage)
            top_number = 0
            for t in range(n):
                top_number += match_set[t][2]
            top_concentration.append([str(i)+'_'+str(j), top_number/len(matches)])
            dz.append(top_number/len(matches))
    
#    fig = plt.figure(figsize= (40, 24))
    fig = plt.figure(figsize= (15, 8))
    ax1 = fig.add_subplot(111, projection='3d')
    
    z3 = np.zeros(len(x3))    
    dx = np.ones(len(x3))*0.5
    dy = np.ones(len(x3))*0.5           
    
    
    ax1.bar3d(x3, y3, z3, dx, dy, dz, color= 'c')
    
    
    ax1.set_xlabel('Length of antibody amino acids')
    ax1.set_ylabel('Length of antigen amino acids')
    ax1.set_zlabel('Relative frequency of the top 5% cores')
    
    ax1.set_xticks(ticks = [1,2,3])
    ax1.set_yticks(ticks = [1,2,3])
    # Save
    os.chdir(results_d)
    plt.savefig('Top_concentration.eps')
    plt.show()
    
    os.chdir(results_d)
    with open('top_concentration', 'w') as f:
        json.dump(top_concentration, f)




def Top_aa_distribution(results_d, t_n=6):
    os.chdir(results_d)
    with open('aa_distribution', 'r') as f:
        aa_distribution = json.load(f)
    top_n = {}
    for i in range(1, 4):
        for j in range(1, 4):
            top_n[str(i)+'_'+str(j)] = {}
    
            for key, value in aa_distribution[str(i)+'_'+str(j)].items():
                value.sort(key = lambda x:x[1], reverse = True)
                top_n[str(i)+'_'+str(j)][key] = value[:t_n]
    with open('top_six', 'w') as f:
        json.dump(top_n, f)

################################################################################
'''
Combination:
    This function combines all the above functions 
Input:
    core_d: the directory of the cores
    results_d: the directory gives all the results
'''
def Combination(core_d, results_d):
    core_distribution = Distribution_over_CDR(core_d)                               
    os.chdir(results_d)            
    with open('core_distribution', 'w') as f:
        json.dump(core_distribution, f)

    distribution_of_aa = Distribution_of_aa()
    os.chdir(results_d)
    with open('aa_distribution', 'w') as f:
        json.dump(distribution_of_aa, f) 
        
    Draw_aa_distribution(results_d)
    
    Top_one_to_one(core_d, results_d)
    
    Cores_distribution_over_CDR(results_d)
    
    Core_distribution_over_match_type(results_d)
    
    Draw_binary_test_AUC(results_d)
    
    Draw_discrimination_AUC(results_d)
    
    Fomalize_results(results_d)
    
    Affinity_pred(results_d)
    
    Top_aa_distribution(results_d, t_n=6)
    
    Top_concentration(core_d, results_d, percentage=0.05)

##########################################################
if __name__ == '__main__':
    core_d = '/home/leo/Documents/Database/Pipeline_New/Complexes/Cores'
    results_d = '/home/leo/Documents/Database/Pipeline_New/Complexes/Results'
    Combination(core_d, results_d)
