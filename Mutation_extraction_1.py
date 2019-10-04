
import numpy as np
import json
import os
import pandas as pd
from pyexcel_ods import get_data, save_data
#import math
os.chdir('/home/leo/Documents/Database/Pipeline_New/Codes')
from AAC_2 import  Get_contact
from FrameConstraint import Get_consecutive
os.chdir('/home/leo/Documents/Database/Pipeline_New/Codes')
from RBFN_coverage import Distance_matrix, Design_matrix

os.chdir('/home/leo/Documents/Database/Data_Code_Publish/Codes')
from Predict_affinity import Predict_affinity, Analyze_resutls

##########################################################

##################################################################################
'''
Coordinates:
    a function to extract coordinates. This one is designed for Affinity validation, and it is talored from the
    'Coordinates' in AAC_2_2 but more efficient.
Input:
    file, combined_chain_id: They are the same as in AAC_2
    mutations_pos: a list gives the mutation position
    mutaion_chain: a string gives the name of the mutation chain
'''
def Coordinates(file, combined_chain_id, mutations_pos, mutation_chain):
    # creat an empty dictionary to contain the results
    cdn = {}
    for i in combined_chain_id[0]:
        cdn['h1'+i], cdn['h2'+i], cdn['h3'+i] = [], [], []
    for i in combined_chain_id[1]:
        cdn['l1'+i], cdn['l2'+i], cdn['l3'+i] = [], [], []
    for i in combined_chain_id[2]:
        cdn[i] = []
        
    # creat a tracker dictionary, and a counter dictionary
    tracker = {}
    counter = {}
    ids = combined_chain_id[0] + combined_chain_id[1] + combined_chain_id[2]
    for i in ids:
        tracker[i] = ''
        counter[i] = -1
        
    # creat a dictionary to indicate the types of chains
    chain_type = {}
    for i in combined_chain_id[0]:
        chain_type[i] = 'H'
    for i in combined_chain_id[1]:
        chain_type[i] = 'L'
    for i in combined_chain_id[2]:
        chain_type[i] = 'A'
        
#    l_range = list(range(1000))
    l_range = []
    l_range.extend(list(range(23, 41)))
    l_range.extend(list(range(49, 64)))
    l_range.extend(list(range(89, 111)))
#        l_range = [[23, 40], [49, 63], [89, 110]]
#        h_range = [[25, 37], [50, 71], [99, 129]]
#    h_range = list(range(1000))
    h_range = []
    h_range.extend(list(range(25, 38)))
    h_range.extend(list(range(50, 72)))
    h_range.extend(list(range(99,130)))
    '''Set the l_range, h_range, and a_range'''
    if mutation_chain in combined_chain_id[0]:
        for pos in mutations_pos:
            if pos not in h_range:
                return []
        h_range = mutations_pos
        l_range = []
        a_range = list(range(10000))
    elif mutation_chain in combined_chain_id[1]:
        for pos in mutations_pos:
            if pos not in l_range:
                return []
        h_range = []
        l_range = mutations_pos
        a_range = list(range(10000))
    elif mutation_chain in combined_chain_id[2]:
        a_range = mutations_pos

    
    # extract the coordinates
    for line in file:
        if line[:6] == 'ATOM  ' and line[21] in ids:
            # update the parameters
            if tracker[line[21]]!= line[22:27]:
                counter[line[21]] += 1
                tracker[line[21]] = line[22:27]
            # extract all the parameters corresponding to line[21]
            c_type = chain_type[line[21]]
            count = counter[line[21]]
            # collect the coordinates according to c_type
            if c_type == 'H':
                if count in h_range:
                    cdn['h1'+line[21]].append([float(line[30:38]), float(line[38:46]), float(line[46:54]),
                                       count, line[17:20]])
    
            if c_type == 'L':
                if count in l_range:
                    cdn['l1'+line[21]].append([float(line[30:38]), float(line[38:46]), float(line[46:54]),
                                       count, line[17:20]])

            if c_type == 'A':
                if count in a_range:
                    cdn[line[21]].append([float(line[30:38]), float(line[38:46]), float(line[46:54]),
                                           count, line[17:20]])            
    
    return cdn


####################################################################
def Sub_get_workable(contact, mutation_match_parameter, pdbid_sequence):
    workable_all = []
    
#    pdbid = mutation_match_parameter['pdbid']
    chain = mutation_match_parameter['mut_chain']
    mutations_pos = mutation_match_parameter['mut_pos']
#    matched_ids = mutation_match_parameter['matched_ids']
#    combined_ids = mutation_match_parameter['combined_ids']
    opposite_chain = mutation_match_parameter['opposite_chain']
    Ab_Ag = mutation_match_parameter['Ab_Ag']
    mut_aa = mutation_match_parameter['mut_aa']
    DDG = mutation_match_parameter['DDG']
    mut_id = mutation_match_parameter['mut_id']
    
    if Ab_Ag == 'Ab':
        mut_chain_pos = 2; opp_chain_pos = 3; aa_pos = 1; opp_aa_pos = 2
    elif Ab_Ag == 'Ag':
        mut_chain_pos = 3; opp_chain_pos = 2; aa_pos = 2; opp_aa_pos = 1

    for fcdn in contact:
        if chain == fcdn[0][mut_chain_pos] and fcdn[0][opp_chain_pos] == opposite_chain:
            for pos, aa in zip(mutations_pos, mut_aa):
                if pos == fcdn[aa_pos]:
                    wt_pair = ['','']
                    wt_pair[aa_pos-1] = [pdbid_sequence[chain][pos]]
                    wt_pair[2-aa_pos] = [pdbid_sequence[opposite_chain][fcdn[opp_aa_pos]]]
                    mut_pair = wt_pair[:]
                    mut_pair[aa_pos-1] = [aa]
                    workable_all.append([wt_pair, mut_pair, fcdn[3], DDG, mut_id]) 
    return workable_all
#####################################################################
'''
Select_contact_opposite:
    a function to select the contact from among all the contact
    and the all the possible opposite
Input:
    mutation_match_parameter:
        a dictionary, contains
        mutation_match_parameter['pdbid']
        mutation_match_parameter['mutation_chain']: 
            a string, gives one of the names of the mutation chains
        mutation_match_parameter['mutations_pos']:
            a list of integers, gives the positions of the mutations on the mutations_chain, it should
            in an increasing order.
        mutation_match_parameter['matched_ids']:
        mutation_match_parameter['combined_ids']:###Can we creat the combined ids according to the matched ids?
            the combined ids corresponding to the pdbid
        mutation_match_parameter['opposite_chain']:
            a string, give one of the names of the opposite chains
        mutation_match_parameter['Ab_Ag']: gives the information of whether the 
                                           mutation chain is an Ab chain or an Ag chain.
       
        cutoff: 
            a float, give the relaxed cutoff value
'''
def One_dict_Workable(mutation_match_parameter, sequence):
    # take out the parameter
    cutoff = mutation_match_parameter['cutoff']
    moving = mutation_match_parameter['moving']
    moving_step = mutation_match_parameter['moving_step']
    moving_start= mutation_match_parameter['moving_start'] 
    moving_end = mutation_match_parameter['moving_end']
    form = mutation_match_parameter['form']
    
    pdbid = mutation_match_parameter['pdbid']
    chain = mutation_match_parameter['mut_chain']
    mutations_pos = mutation_match_parameter['mut_pos']
    matched_ids = mutation_match_parameter['matched_ids']
    combined_ids = mutation_match_parameter['combined_ids']

    
    # Extract the required data
    '''************************************'''
    os.chdir('/home/leo/Documents/Database/Data_Code_Publish/Structures/imgt')
    '''*******************************************************'''
    
    with open(pdbid+'.pdb', 'r') as file:
        cdn = Coordinates(file, combined_ids,  mutations_pos, chain)
    # Check if the cdn is []
    if cdn == []:
        return []
#        pdbid_sequence = Chain_seq(file, combined_ids)
    pdbid_sequence = sequence[pdbid]
    # store the paires in paires
    workable_all = []
    if not moving:    
        contact = Get_contact(cdn, matched_ids, cutoff)       
        if contact != []:
            workable_all = Sub_get_workable(contact, mutation_match_parameter, pdbid_sequence)
                         
                   
    if moving:
        moving_cutoff = moving_start
        while workable_all == [] and moving_cutoff <= moving_end:           
            contact = Get_contact(cdn, matched_ids, moving_cutoff)       
            if contact != []:
                workable_all = Sub_get_workable(contact, mutation_match_parameter, pdbid_sequence)
         
            moving_cutoff += moving_step
            
    # Truncate according to the form       
    if workable_all != []:
        if form == 'one':
            workable_all.sort(key=lambda x:x[2], reverse=True)
            workable_one_dict = [workable_all[0]]
        elif form == 'multiple':
            workable_one_dict = workable_all[:]
    else:
        workable_one_dict = workable_all
        
    return workable_one_dict


##################################################################
def Parse_mutations(mut_data, matched_ids):
    parsed_list = []
    for pdbid, values in mut_data.items():
        matched_ids_pdb = matched_ids[pdbid]
        # GET THE DATA FOR ONE MUTATION SET
        for i, mu in enumerate(values):
            if mu[0] == pdbid:
                start = i
                DDG = 1000
                mut_id = mu[5]
            elif mu[0] == 'affinities':
                end = i
                DDG = mu[2]
            ###########################################
            if DDG != 1000:
                mut_chains = []; one_parsed_list = []
                for mut in values[start:end]:
                    if mut[1] not in mut_chains:
                        mut_chains.append(mut[1])
                for chain in mut_chains:
                    # find the opposite chains
                    for matched in matched_ids_pdb:
                        for k, v in enumerate(matched):
                            if chain in v:
                                if k == 0 or k == 1:
                                    Ab_Ag = 'Ab'
                                    opposite_chains = matched[2]
                                else:
                                    Ab_Ag = 'Ag'
                                    opposite_chains = matched[0]+matched[1]
                    # find the mut aa and mut pos
                    mut_aa = []; mut_pos = []
                    for mut in values[start:end]:
                        if mut[1] == chain:
                            mut_aa.append(mut[3])
                            mut_pos.append(mut[2])

                    one_parsed_list.append([pdbid, Ab_Ag, chain, mut_pos, mut_aa, opposite_chains, DDG, mut_id])
                parsed_list.append(one_parsed_list)
                    
    return parsed_list
###############################################################################
def Get_workable(matched_ids,combined_ids,sequence, parsed_list, cutoff, moving,\
                 moving_step, moving_start, moving_end, form):
    workable = []
    mutation_match_parameter = {}
    mutation_match_parameter['cutoff'] = cutoff 
    mutation_match_parameter['moving'] = moving 
    mutation_match_parameter['moving_step'] = moving_step 
    mutation_match_parameter['moving_start'] = moving_start
    mutation_match_parameter['moving_end'] = moving_end 
    mutation_match_parameter['form'] = form 
    for one_parsed_list in parsed_list:
        one_mut_workable = []
        if one_parsed_list != []:
            for one_chain in one_parsed_list:
                print('Working on   ', one_chain[7])
                # Load the mutaion_mut_parameter
                mutation_match_parameter['pdbid'] = one_chain[0]
                mutation_match_parameter['Ab_Ag'] = one_chain[1]
                mutation_match_parameter['mut_chain'] = one_chain[2]
                mutation_match_parameter['mut_pos'] = one_chain[3]
                mutation_match_parameter['mut_aa'] = one_chain[4]
                mutation_match_parameter['DDG'] = one_chain[6]
                mutation_match_parameter['mut_id'] = one_chain[7]
                mutation_match_parameter['matched_ids'] = matched_ids[one_chain[0]]
                mutation_match_parameter['combined_ids'] = combined_ids[one_chain[0]]
                for opp_chain in one_chain[5]:
                    mutation_match_parameter['opposite_chain'] = opp_chain              
                    one_chain_workable = One_dict_Workable(mutation_match_parameter, sequence)
                    one_mut_workable.extend(one_chain_workable)
                    
        if one_mut_workable != []:      
            workable.append(one_mut_workable)
    
    return workable
###############################################################################
if __name__ == '__main__':
    # Read data
    os.chdir('/home/leo/Documents/Database/Data_Code_Publish/Structures')
    with open('combined_ids', 'r') as f:
        combined_ids = json.load(f)
    with open('matched_ids', 'r') as f:
        matched_ids = json.load(f)
    with open('sequence', 'r') as f:
        sequence = json.load(f)
    # parse the mutaion
    os.chdir('/home/leo/Documents/Database/Pipeline_New/Codes/Results')
    mut_data = get_data('Keating.ods')
    
    parsed_list = Parse_mutations(mut_data, matched_ids)
    
    # Get the workable
    cutoff = 6
    moving = True
    moving_step = 0.25
    moving_start = 3
    moving_end= 8
    form = 'one'
    workable = Get_workable(matched_ids,combined_ids,sequence, parsed_list, cutoff, moving,\
                 moving_step, moving_start, moving_end, form)



working_d = '/home/leo/Documents/Database/Data_Code_Publish/Codes/Results'
pred_res_old = Predict_affinity(workable, working_d, binary = True)
selected_cut_DDG, AUC, TPR, FPR, correct_ratio = \
            Analyze_resutls(pred_res_old, cut_DDG_lower=0.5, cut_DDG_upper=100)



#os.chdir('/home/leo/Documents/Database/Pipeline_New/Complexes/Results')
#with open('one_to_one_affinity', 'r') as f:
#    one_to_one_affinity = json.load(f)
#
#one_to_one_affinity['affinity_results']['results_dict'].keys()    
#one_to_one_affinity.keys()
#one_to_one_affinity['0.5_1000']['auc_all']
## Extract the information
#results = one_to_one_affinity['affinity_results']['results_dict']['single_True_mode_single_moving_True_binary_True']['0']
#workable_old = []
#mut_sets=[]; DDG = []; mut_id = []
#for res in results:
#    one_mut_set = []; one_workable = []
#    for i in range(1, len(res)-1):
#        if res[i] != []:
#            one_mut_set.append([res[i][0][0], res[i][1][0]]) 
#    if one_mut_set != []:
#        for mut in one_mut_set:
#            one_workable.append([mut[0],mut[1], 1, res[0][2][2], res[0][3]]) 
#        mut_sets.append(one_mut_set)
#        DDG.append(res[0][2][2])
#        mut_id.append(res[0][3])
#        workable_old.append(one_workable)


