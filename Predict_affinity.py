'''PREDICT MUTATION'''
import json
import xlrd
import os
import copy
from pyexcel_ods import get_data, save_data
from Extract_mutation_pairs import Formalized_contacting
'''#############################################################'''
# First change the mutation into the inputs of function Complete_extract_contacting
# The basic strategy is to process each mutaion set with measured affinity as a unit
# and iterate all the mutaion set
'''  
search_para = {}          
search_para['moving'] 
search_para['step_size']
search_para['start_dist']
search_para['end_dist']
search_para['cut_dist']
search_para['form']
search_para['within_range']
search_para['pdbid']
earch_para['mut_pos']
search_para['mut_chain_id']
'''

'''##############################################################
Processing_mutation_set: a sub function of Workable input
'''

def Processing_mutation_set(mutation_set):
    # Separate the information
    DDG = mutation_set[-1][1]
    pdbid = mutation_set[0][0]
    mut_id = mutation_set[0][5]
    mutations = mutation_set[:-1]
    
    # Get the chains, and pos
    chains = []
    for mu in mutations:
        if mu[1] not in chains:
            chains.append(mu[1])
    
    mut_info = []
    for chain in chains:
        pos = [mu[2] for mu in mutations if mu[1] == chain]
        mut_aa = [mu[3] for mu in mutations if mu[1] == chain]

        mut_info.append(copy.deepcopy([pdbid, chain, pos, mut_aa, DDG, mut_id]))
    
    return mut_info
# Read the parsed mutation information
'''
Output:
    workable_input:
        a list with elements in the form 
        
        [...,
        [pdbid, mut_chain_id, [pos], [mut_aa], affinities, mut_id],
        ...]
        
        The element is a list with the number of elements equals the number of 
        mutation chains.
'''
def Workable_input(mutation_d):
    os.chdir(mutation_d)
    formated = get_data('Formated.ods')

    workable_input = []
    for pdb, mutations in formated.items():
        begin = 0; end = 0
        for ind, mut in enumerate(mutations):
            if mut[0] == 'affinities':
                end = ind + 1
                mutation_set = mutations[begin: end]
                begin = end
                # processing the mutation set
                mut_info = Processing_mutation_set(mutation_set)
                workable_input.append(copy.deepcopy(mut_info))
            
    return workable_input
        
mutation_d = '/home/leo/Documents/Database/Data_Code_Publish/Mutations'
workable_input = Workable_input(mutation_d)   
'''##############################################################3'''
'''
Generate_wt_mut_one: a function to generat wt_mut pairs for one_tuple
Input:
    one_tuple:
        in the form
        
        ( mut_pos, mut_aa, mut_chain_id,  op_pos, op_aa, op_chain_id,\
        contact_number, mut_chain_type, op_chain_type)
        
        if form == 'flanked', mut_aa and op_aa are lists
Output:
    Ab_Ag_wt_mut_pair: a list in the form
    [[wt_Ab_Ag_pair], [mut_Ab_Ag_pair]]
'''
def Generate_wt_mut_one(one_tuple, form, mut_aa):
    mut_chain_type = one_tuple[7]
    # Extract wt information from one_tuple and generate mut_info
    if form == 'one' or 'multiple':
        wt_pair = [[one_tuple[1]], [one_tuple[4]]]
        mut_pair = [[mut_aa], [one_tuple[4]]]

    elif form == 'flanked':
        wt_pair = [one_tuple[1], one_tuple[4]]
        
        tri = one_tuple[1]
        tri[1] = mut_aa
        mut_pair = [tri, one_tuple[4]]
        
    # Put the Ab first
    if mut_chain_type == 'A':
        Ab_Ag_wt_pair = [wt_pair[1], wt_pair[0]]
        Ab_Ag_mut_pair = [mut_pair[1], mut_pair[0]]
    
    Ab_Ag_wt_mut_pair = [Ab_Ag_wt_pair, Ab_Ag_mut_pair]
    
    return Ab_Ag_wt_mut_pair
'''
One_chian_output: a sub fuction of Workable_output
Output:
    a list with elements in the following form: 
    
    [[wt_Ab_Ag_pair], [mut_Ab_Ag_pair], DDG, mut_id]
    
'''
def One_chian_output(one_chain, search_para, combined_ids, sequence, structure_d):
    
    search_para['pdbid'] = one_chain[0]
    search_para['mut_chain_id'] = one_chain[1]
    search_para['mut_pos'] = one_chain[2]
    
    mut_pos_list = one_chain[2]
    mut_aa_list = one_chain[3]
    DDG = one_chain[4]
    mut_id = one_chain[5]
    
    form = search_para['form']
    
    formalized_pairs = \
        Formalized_contacting(search_para, combined_ids, sequence, structure_d)
    
    # reformat the formalized pairs with more information and transform into wt mut pairs
    Ab_Ag_wt_mut_pairs_list = []
    for one_tuple in formalized_pairs:
        for pos, aa in zip(mut_pos_list, mut_aa_list):
            if one_tuple[0] == pos:
                mut_aa = aa
                Ab_Ag_wt_mut_pair = Generate_wt_mut_one(one_tuple, form, mut_aa)
                # Attach the DDG and the mutid
                Ab_Ag_wt_mut_pair.extend([DDG, mut_id])
                # Load to the Ab_Ag_wt_mut_pairs_list
                Ab_Ag_wt_mut_pairs_list.append(copy.deepcopy(Ab_Ag_wt_mut_pair))

        
    return Ab_Ag_wt_mut_pairs_list
'''
Input:
    search parameter contains the following values:
        search_para['moving'] 
        search_para['step_size']
        search_para['start_dist']
        search_para['end_dist']
        search_para['cut_dist']
        search_para['form']
        search_para['within_range']
        
Output:
    workable, a list with elements lists in the following form
    
    [...'
    [[wt_Ab_Ag_pair], [mut_Ab_Ag_pair], DDG, mut_id],
    ...]
    
    The length of this list is the same as the number of muation chains is a mutation_set:
        mutaion_set is defined as a set of mutations of which the affinity of the complex
        is measured.
'''
def Workable_output(workable_input, search_para, combined_ids, sequence, structure_d):
    workable = []

    for one_mutation in workable_input:
        one_mut_workable = []
        for one_chain in one_mutation:
            Ab_Ag_wt_mut_pairs_list =\
                One_chian_output(one_chain, search_para, combined_ids, sequence, structure_d)
            one_mut_workable.extend(Ab_Ag_wt_mut_pairs_list)
            
        workable.append(copy.deepcopy(one_mut_workable))
    
    return workable

'''############################################################################'''
def Predict_affinity():
    pass
def Analyze_resutls():
    pass







































