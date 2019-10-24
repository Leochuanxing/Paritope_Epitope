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




















