'''
This file is to select the complexes with given resolution
'''
import os
import json
import copy
import math
import pandas as pd
import numpy as np
'''########################################################################################'''
class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)
'''########################################################################################'''        
def Number_date(date):
    date_number = ''
    date_list = date.split('/')
    if int(date_list[2]) < 20:
        date_number += '20'
    else: date_number += '19'
    date_number += date_list[2]
    date_number += date_list[0]
    date_number += date_list[1]
        
    date_number = int(date_number)
    
    return date_number
'''#################################################################################'''
def Get_matched_ids(summary):
    matched_ids = {}
    for i in range(summary.shape[0]):
        pdb = summary.iloc[i]['pdb']
        if summary.iloc[i]['Hchain'] != summary.iloc[i]['antigen_chain'] and \
            summary.iloc[i]['Lchain'] != summary.iloc[i]['antigen_chain']:
            if pdb not in matched_ids:
                matched_ids[pdb] = [[summary.iloc[i]['Hchain'], summary.iloc[i]['Lchain'],\
                            summary.iloc[i]['antigen_chain']]]
            else:
                matched_ids[pdb].append([summary.iloc[i]['Hchain'], summary.iloc[i]['Lchain'],\
                            summary.iloc[i]['antigen_chain']])
    return matched_ids
'''##################################################################################'''
'''
Id_dict, is to extract the names of the heavy chain, light chain and the antigen
        chain for each pdb id from the "summary" file.
Input: 
    file, is the summary file
Output:
    id_dict, a dictionary with the keys pdb ids and the values in the form of 
        lists. Take the complex "4ydl" for example, the output is 
        [['B', 'C', 'A'], ['H', 'L', 'G']], where ['B', 'C', 'A'] 
        means the heavy chain is B, the light chain is C and the antigen chain 
        is A. And ['H', 'L', 'G'] has the similar meanings. 
        From this result, we can tell that this complex may be composed of two 
        subunits, and those two units may be the same.    
'''
      
def Id_dict (file):
    id_dict = {}
    for l in file:
#        Make sure the line is long enough
        if len(l) >= 16:
            a = l.split('\t')
#        Deal with the | in a[4]            
            for i in a[4].split('|') :
                temp = [a[1].strip(), a[2].strip(), i.strip()]
                for i in range(0,3):
                    if temp[i] == 'NA' or not temp[i].isupper():
                        temp[i] =''                       
                if a[0].strip() in id_dict and temp[2] != '':
                    id_dict[a[0].strip()].append (temp)
                elif a[0].strip() not in id_dict and temp[2] != '':
                    id_dict[a[0].strip()] = [temp]                    
    return id_dict 
'''###################################################################################'''         
'''
Combined_chian_id_dict is to combined the ids of the heavy chain, light chain and 
        antigen chain together, to facilitate the iteration.
Input: 
    id_dict, which is the output of the function Id_dict above
Output:
    combined_chain_id_dict, the output will be in the form of a list.
    Take the complex "4ydl" for example, the output is 
    ['BH', 'CL', 'AG'], "B" and "H" are heavy chains, "C" and "L" are light chains,
    "A" and "G" are antigen chains
'''
def Combined_chain_id_dict (id_dict):
    combined_chain_id_dict = {}
    for i in id_dict:
        temp = ['' ,'' ,'' ]        
        for j in id_dict[i]:
            temp = [temp[0]+j[0], temp[1]+j[1], temp[2]+j[2]]
        combined_chain_id_dict[i] = temp
    return combined_chain_id_dict   

########################################################################
#Eliminate the wierd pdb, and set aside those wierd pdb ids
# by weird, it means a chain id shows up more than once
'''
witch_hunt,. is to find out the pathological complex with the same chain shows up
          both as the antigen chain and the antibody chain.
inputs: conbined_chain_id, a dictionary with keys pdb id, and values in the form
        of ['BDF', 'ACE', 'GH']
return: witches, a list of weird pdb ids
        GoodPeople, a dictionary with keys of good pdb ids, and values in the form
        of ['BDF', 'ACE', 'GH']
'''
def Eliminate_abnomal(conbined_chain_id, id_dict):
    # creat a container to contain the witches and GoodPeople
    # the witches are the ones with the same chain appears both as the antigen chain
    # and as the antibody chain.
    witches = []
    GoodPeople = conbined_chain_id
    for pdbid in conbined_chain_id:
        antibody_ids = conbined_chain_id[pdbid][0] + conbined_chain_id[pdbid][1] 
        antigen_ids = conbined_chain_id[pdbid][2]
        witch = False
        # here is to find the witch
        for chain in antibody_ids:
            if chain in antigen_ids:
                witch = True
                break
        if witch:
            witches.append(pdbid)
    for witch in witches:
        del GoodPeople[witch]
        del id_dict[witch]
        
            
    return witches, GoodPeople, id_dict
# Hunt the witches and generate the good_combined_ids  and good_matched_ids for 


########################################################################
def Combined_Matched_ids(working_directory):
    os.chdir(working_directory)
    # Work on the summary of complexes fo peptide antigen
    summary = open('summary_peptide_antigen_4A.tsv', 'r')# case sensitive
    files = summary.readlines()        
    summary.close

    id_dict = Id_dict (files)
    combined_chain_id_dict = Combined_chain_id_dict (id_dict)
    witches, good_combined_ids, good_matched_ids = Eliminate_abnomal(combined_chain_id_dict, id_dict)
    
    # Work on the summary of complexes of protein antigen
    summary = open('summary_protein_antigen_4A.tsv', 'r')# case sensitive
    files = summary.readlines()        
    summary.close

    id_dict = Id_dict (files)
    combined_chain_id_dict = Combined_chain_id_dict (id_dict)
    witches, good_combined_ids2, good_matched_ids2 = Eliminate_abnomal(combined_chain_id_dict, id_dict)
    
    # Combine the two dictionaries
    good_combined_ids.update(good_combined_ids2)
    good_matched_ids.update(good_matched_ids2)
    
    return good_combined_ids, good_matched_ids

#working_directory = '/home/leo/Documents/Database/Data_Code_Publish/Structures'
#good_combined_ids, good_matched_ids = Combined_Matched_ids(working_directory)

'''#####################################################################################3'''
'''

Generate_testing_training_ids:
    Generate the testing_combined_ids, testing_matched_ids, training_combined_ids, 
    training_matched ids
Inputs:
    cut_resolution: only select the complexes with the resolution no larger than cut_resolution
    precentage: the percentage of complexes will be used as the testing set
Outputs:


'''
def Generate_testing_training_ids(cut_resolution=3, precentage=0.1):
   
    os.chdir('/home/leo/Documents/Database/Data_Code_Publish/Structures')
    summary1 = pd.read_table('summary_peptide_antigen_4A.tsv', sep='\t')
    summary2 = pd.read_table('summary_protein_antigen_4A.tsv', sep='\t')
    
    summary = pd.concat([summary1, summary2], axis = 0, ignore_index=True)


    # Arrange the complex according to the order of submission
    number_date = summary['date'].apply(Number_date)
    # Attach this number_date series to the dataframe
    summary['number_date'] = number_date
    # Arrange the rows according to the number_data
    summary_sorted = summary.sort_values('number_date')

    # Select the cut_resolution
    summary_sorted_cut = summary_sorted[summary_sorted['resolution']<= cut_resolution]

    # Generate the combined ids and the matched ids
    pdbids_date = []
    for i in range(summary_sorted_cut.shape[0]):
        pdb_d = [summary_sorted_cut.iloc[i]['pdb'], summary_sorted_cut.iloc[i]['number_date']]
        if pdb_d not in pdbids_date:
            pdbids_date.append(pdb_d)        

    # Get the training pdbids and the testing pdbids
    percentage = 0.1
    n_train = math.floor((1-percentage) * len(pdbids_date))
    training_ids = pdbids_date[:n_train]
    testing_ids = pdbids_date[n_train:]
    cut_date = pdbids_date[n_train][1]
#        # Get the cut date
#        training_summary = summary_sorted_cut[summary_sorted_cut['number_date'] <= cut_date]
#        testing_summary = summary_sorted_cut[summary_sorted_cut['number_date'] > cut_date]
# Get the combined ids and the matched ids
    return training_ids, testing_ids, cut_date

'''##########################################################################################'''
def main(working_directory, cut_resolution=3, precentage=0.1):

    working_directory = '/home/leo/Documents/Database/Data_Code_Publish/Structures'
    good_combined_ids, good_matched_ids = Combined_Matched_ids(working_directory)
    
    training_ids_date, testing_ids_date, cut_date = Generate_testing_training_ids(cut_resolution=3, precentage=0.1)
    
    # Load the ids and dates to a dictionary
    train_test_ids_dates = {}
    train_test_ids_dates['training_ids_dates'] =[[x[0], int(x[1])] for x in training_ids_date]
    train_test_ids_dates['testing_ids_dates'] = [[x[0], int(x[1])] for x in testing_ids_date]

    
    
    # Load the selected 
    training_combined_ids={}; training_matched_ids = {}
    testing_combined_ids = {}; testing_matched_ids = {}
    for train_pdb_date in training_ids_date: 
         if train_pdb_date[0] in good_matched_ids:
            training_combined_ids[train_pdb_date[0]] = good_combined_ids[train_pdb_date[0]]
            training_matched_ids[train_pdb_date[0]] = good_matched_ids[train_pdb_date[0]]
        
    for test_pdb_date in testing_ids_date:
        if test_pdb_date[0] in good_matched_ids:
            testing_combined_ids[test_pdb_date[0]] = good_combined_ids[test_pdb_date[0]]
            testing_matched_ids[test_pdb_date[0]] = good_matched_ids[test_pdb_date[0]]

    # Save the results
    with open('matched_ids', 'w') as f:
        json.dump(good_matched_ids, f)
    with open('combined_ids', 'w') as f:
        json.dump(good_combined_ids, f)
    with open('training_combined_ids', 'w') as f:
        json.dump(training_combined_ids, f)
    with open('training_matched_ids', 'w') as f:
        json.dump(training_matched_ids, f)
    with open('testing_combined_ids', 'w') as f:
        json.dump(testing_combined_ids, f)
    with open('testing_matched_ids', 'w') as f:
        json.dump(testing_matched_ids, f)
        
    with open('train_test_ids_dates', 'w') as f:
        json.dump(train_test_ids_dates, f)
##        
    return 

if __name__ == '__main__':
    working_directory = '/home/leo/Documents/Database/Data_Code_Publish/Structures'
    main(working_directory, cut_resolution=3, precentage=0.1)
