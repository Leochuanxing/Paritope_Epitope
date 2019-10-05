import os
import json
import copy
import sys
from pyexcel_ods import get_data
'''################################################################'''
'''
Sub_extract_cdn: A function to extract the coordinates for given parameters
Input: 
    file: relevant pdb_file
    chain_id: str with one letter
    pos: a list os relevant indices
Output:
    cdn: a list. Each element in the list is in the form of 
    [pos, aa, x, y, z, chain_id]
'''
def Sub_extract_cdn(file, chain_id, pos):
    
    cdn = []    
    # creat a tracker and a counter
    tracker = ''; counter = -1
    
    for line in file:
        if line[:6] == 'ATOM  ' and line[21] == chain_id:
            # update the parameters
            if tracker!= line[22:27]:
                counter += 1; tracker = line[22:27]
            if counter in pos:
                cdn.append([counter, line[17:20], float(line[30:38]), float(line[38:46]), float(line[46:54]), chain_id])
                    
    return cdn
'''##############################################################'''
'''
Inputs:
    mut_pos: a list, gives the mutation positions
    mut_chain_id: a one letter string, gives the name of the chain where the mutation is located
'''
def Extract_cdn(pdbid, mut_pos, mut_chain_id, combined_ids, sequence, structure_d):
    # Decide the mut_chain_type      
    combined = combined_ids[pdbid]
#    for matched in combined_pdb:
    if mut_chain_id in combined[0]:
        mut_chain_type = 'H'
        op_chain_ids = combined[2]
    elif mut_chain_id in combined[1]:
        mut_chain_type = 'L'
        op_chain_ids = combined[2]
    elif mut_chain_id in combined[2]:
        mut_chain_type = 'A'
        op_chain_ids = combined[0]+combined[1]
        
        
    # Extract the coordinates for given pdbid, chain_id, pos
    os.chdir(structure_d)
    # Extract the cdn ofthe mutation
    with open(pdbid + '.pdb', 'r') as file:
        mut_cdn = Sub_extract_cdn(file, mut_chain_id, mut_pos)
    file.close
    # Extract the cdn of the opposite amino acids
    op_pos = list(range(10000))
    op_cdn =[]
    for op_chain in op_chain_ids:
        with open(pdbid+'.pdb', 'r') as file:
            op_cdn.extend(Sub_extract_cdn(file, op_chain, op_pos))
            file.close
            
    return mut_cdn, mut_chain_type, op_cdn, op_chain_ids
'''#####################################################################'''
'''
Output:
    contacting_pairs: a list, with the elements in the folowing form
    [ mut_pos, mut_aa, mut_chain_id,  op_pos, op_aa, op_chain_id, contact_number, mutation_chain_type]
    for example:
    [ 2, TYR, H, 8, SER, A,  17, 'H']
'''
def Sub_extract_contacting_pairs(mut_cdn, op_cdn, mut_chain_type, dist):
    
    contact_pairs = []
    
    contact = []
    dist_square = dist * dist
    for mut in mut_cdn:
        for op in op_cdn:
            d_square = 0
            for x, y in zip(mut[2:5], op[2:5]):
                d_square += (x-y) * (x-y)
                if d_square > dist_square:
                    break
            if d_square <= dist_square:
                contact.append((mut[0], mut[1], mut[5], op[0], op[1], op[5]))
                
    if contact != []:
        # Separate all the coordinates by chains
        for six in set(contact):
            cn = 0
            for s in contact:
                if s == six:
                    cn += 1
                    
            # Extend the tuple        
#            six_c = copy.deepcopy(six)
            six_c = list(six)
            six_c.extend([cn, mut_chain_type])
            six_c = tuple(six_c)
            # load
            contact_pairs.append(six_c)
            
    return contact_pairs
    
'''#####################################################################'''
'''
Inputs:
    search_para: a dictionary, contains
        ['moving'] = True or False
        ['step_size']: a float, defines the step size for each moving step
        ['start_dist']: a float, defines the starting distance of the moving method
        ['end_dist']: a float, defines the ending distance of the moving method
        ['cut_dist']: a float, defines the cut off distance
        ['form']: string, takes the value in ['one', 'multiple', 'flanked']
        ['within_range']: boolean, True or False. If True, means the Ab amino acide in the contacting pair 
                         should be within the range fo the cdrs
        
        Ff ['moving'] is True, ['step_size'], ['start_dis'] and ['end_dist'] have to be 
        assigned. If['moving'] is False, ['cut_dist'] has to be defined.
        
        For convenience, we list the attributes bellow:
            search_para['moving'] 
            search_para['step_size']
            search_para['start_dist']
            search_para['end_dist']
            search_para['cut_dist']
            search_para['form']
            search_para['within_range']
            
    mut_cdn: The return value of Extract_cdn
    mut_chain_type: The return value of Extract_cdn
    op_cdn: The return value of Extract_cdn
Output:
    contact_pairs: it is a concept related to the mutations introduced in one chain

'''
def Extract_contacting_pairs(mut_cdn, mut_chain_type, op_cdn, search_para):
    # Take out the values of the parameters
    moving = search_para['moving']
    if moving:
        step_size = search_para['step_size']
        start_dist = search_para['start_dist']
        end_dist = search_para['end_dist']
        
        dist = start_dist
        while dist <= end_dist:
            contact_pairs = Sub_extract_contacting_pairs(mut_cdn, op_cdn, mut_chain_type, dist)
            dist += step_size
            if contact_pairs != []:
                break
            
    if not moving:
        cut_dist = search_para['cut_dist']
        contact_pairs = Sub_extract_contacting_pairs(mut_cdn, op_cdn, mut_chain_type, cut_dist)
    
    return contact_pairs

'''#########################################################################'''
'''
search_para:
    ['form']: string, takes the value in ['one', 'multiple', 'flanked']
    
        'one': in a contact_pairs only one 1_to_1 interaction with the largest contact 
                number is selected
        'multiple': in a contact_pairs, all the 1_to_1 interaction are selected and each one
                will be weighted by the contact number
        'flanked': in a contact_pairs, the 1_to_1 pair with the largest contact number is selected
                and the amino acid before and after this contacting pairs will be selected. Therefore
                the 1_to_1 contact pair is converted to a 3_to_3 contact pair.
    ['within_range']: boolean, True or False. If True, means the Ab amino acide in the contacting pair 
            should be within the range fo the cdrs
Output:
    formalized_pairs: a list with element in the form:
        
        (mut_pos, mut_aa, mut_chain_id,  op_pos, op_aa, op_chain_id, contact_number, mut_chain_type, op_chain_type)
        
        For example:
            
        [(31, 'TYR', 'W', 107, 'TRP', 'H', 3, 'A', 'H'),
         (31, 'TYR', 'W', 105, 'SER', 'H', 10, 'A', 'H')]
    
'''
def Formalize_contact_pairs(pdbid, sequence, contact_pairs, search_para, combined_ids):
    
    # Extract different forms#!/usr/bin/env python3
    form = search_para['form']
    
    try:
        if form == 'one':
            formalized_pairs = Form_one(contact_pairs)
        elif form == 'multiple':
            formalized_pairs = Form_multiple(contact_pairs)
        elif form == 'flanked':
            formalized_pairs = Form_flanked(contact_pairs, sequence, pdbid)
    except UnboundLocalError:
        print('Can not find the input form')
        sys.exit()
    
    # Attach the op_chain_type in case it will be used latter
    formalized_pairs_c = copy.deepcopy(formalized_pairs)
    formalized_pairs = []
    combined = combined_ids[pdbid]
    for pair in formalized_pairs_c:
        if pair[5] in combined[0]:
            op_chain_type = 'H'
        elif pair[5] in combined[1]:
            op_chain_type = 'L'
        elif pair[5] in combined[2]:
            op_chain_type = 'A'
        # Attach
        pair_c = list(pair)
        pair_c.append(op_chain_type)
        formalized_pairs.append(tuple(pair_c))
        

    # Check if it is within range, use the same range as in AAC_2
    l_range = list(range(23, 40))
    l_range.extend(list(range(49, 63)))
    l_range.extend(list(range(89,110)))
#    l_range = [[23, 40], [49, 63], [89, 110]]
    h_range = list(range(25, 37))
    h_range.extend(list(range(50, 71)))
    h_range.extend(list(range(99, 129)))
#    h_range = [[25, 37], [50, 71], [99, 129]]
    within_range = search_para['within_range']
    
    if within_range and formalized_pairs != []: 
        for pair in formalized_pairs:               
            if pair[7] == 'H' and pair[0] not in h_range:
                return []
            if pair[7] == 'L' and pair[0] not in l_range:
                return []
            if pair[8] == 'H' and pair[3] not in h_range:
                return []
            if pair[8] == 'L' and pair[3] not in l_range:
                return []

    return formalized_pairs
'''##########################################################################3'''
def Form_one(contact_pairs):
    form_one = []
    if contact_pairs !=[]:
        contact_pairs.sort(key=lambda x : x[6], reverse = True)
        form_one = [contact_pairs[0]]
    return form_one

def Form_multiple(contact_pairs):
    return contact_pairs

def Form_flanked(contact_pairs, sequence, pdbid):
    flanked = []
    
    form = Form_one(contact_pairs)
#    form = Form_multiple(contact_pairs)
    
    if form != []:
        for i in range(len(form)):
            form_one = form[i]
            
            mut_pos = form_one[0]
            mut_chain_id = form_one[2]
            op_pos = form_one[3]
            op_chain_id = form_one[5]
            
            # Find the amino acids before and after
            mut_before_pos = max(0, mut_pos - 1)
            op_before_pos = max(0, op_pos - 1)
            
            mut_after_pos = min(mut_pos + 2, len(sequence[pdbid][mut_chain_id]))
            op_after_pos = min(op_pos + 2, len(sequence[pdbid][op_chain_id]))
            
            flanked_mut = []
            for m in range(mut_before_pos, mut_after_pos):
                flanked_mut.append(sequence[pdbid][mut_chain_id][m])
            if mut_before_pos > mut_pos - 1:
                flanked_mut.insert(0, '')
            if mut_after_pos < mut_pos + 2:
                flanked_mut.append('')
            
            flanked_op = []
            for o in range(op_before_pos, op_after_pos):
                flanked_op.append(sequence[pdbid][op_chain_id][o])
            if op_before_pos > op_pos - 1:
                flanked_op.insert(0, '')
            if op_after_pos < op_pos + 2:
                flanked_op.append('')
                
            # Replace the mut_aa and op_aa in form one with the flanked
            flanked_one = list(form_one)
            flanked_one[1] = flanked_mut
            flanked_one[4] = flanked_op
            flanked.append(tuple(flanked_one))
#            flanked = [tuple(flanked)]# Change back into the same form
        
    return flanked
'''#######################################################################'''
'''
Formalized_contacting:
    This function is for the convenience of using by other modules
Input:
    combined_ids, sequence: as before
    structure_d: the directory of the structures
    search_para: a dictionary containing the following values
            search_para['moving'] 
            search_para['step_size']
            search_para['start_dist']
            search_para['end_dist']
            search_para['cut_dist']
            search_para['form']
            search_para['within_range']
            search_para['pdbid']
            search_para['mut_pos']
            search_para['mut_chain_id']
Output: 
    formalized_pairs, as explained in 'Formalize_contact_pairs' 

REMARK: 
    1.  search_para['pdbid'], a string with length 1
    2.  earch_para['mut_pos'], a list gives the positions of the mutations
'''

def Formalized_contacting(search_para, combined_ids, sequence, structure_d):
    pdbid = search_para['pdbid']
    mut_pos = search_para['mut_pos']
    mut_chain_id = search_para['mut_chain_id']
    # Extract the cdn
    mut_cdn, mut_chain_type, op_cdn, op_chain_ids = \
        Extract_cdn(pdbid, mut_pos, mut_chain_id, combined_ids, sequence, structure_d)
        
    contact_pairs =\
       Extract_contacting_pairs(mut_cdn, mut_chain_type, op_cdn, search_para)
   
    formalized_pairs =\
       Formalize_contact_pairs(pdbid, sequence, contact_pairs, search_para, combined_ids)
    
    return formalized_pairs
'''#############################################################################'''
'''##############################################################
Processing_mutation_set: a sub function of Workable input
'''

def Processing_mutation_set(mutation_set):
    mut_info = []
    # Separate the information
    DDG = mutation_set[-1][1]
    pdbid = mutation_set[0][0]
    mut_id = mutation_set[0][5]
    mutations = mutation_set[:-1]
    
    for mu in mutations:
        mut_info.append([pdbid, mu[1], [mu[2]], [mu[3]], DDG, mut_id])
    
    return mut_info
'''#############################################################################'''
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
                workable_input.append(mut_info)
            
    return workable_input
  
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
    if form == 'one' or form == 'multiple':
        wt_pair = [[one_tuple[1]], [one_tuple[4]]]
        mut_pair = [[mut_aa], [one_tuple[4]]]

    elif form == 'flanked':
        wt_pair_1 = [x for x in one_tuple[1] if x != '']
        wt_pair_2 = [x for x in one_tuple[4] if x != '']
        wt_pair = [wt_pair_1 , wt_pair_2]
        
        tri = one_tuple[1][:]
        tri[1] = mut_aa
        mut_pair_1 = [x for x in tri if x != '']
        mut_pair_2 = [x for x in one_tuple[4] if x!= '']
        mut_pair = [mut_pair_1 , mut_pair_2]
        
    # Put the Ab first
    if mut_chain_type == 'A':
        Ab_Ag_wt_pair = [wt_pair[1], wt_pair[0]]
        Ab_Ag_mut_pair = [mut_pair[1], mut_pair[0]]
    else:
        Ab_Ag_wt_pair = wt_pair
        Ab_Ag_mut_pair = mut_pair
        
    
    Ab_Ag_wt_mut_pair = [Ab_Ag_wt_pair, Ab_Ag_mut_pair]
    
    return Ab_Ag_wt_mut_pair
'''###################################################################################'''

'''
One_chian_output: a sub fuction of Workable_output
Output:
    a list with elements in the following form: 
    
    [[wt_Ab_Ag_pair], [mut_Ab_Ag_pair],contact_number, DDG, mut_id]    
'''
def Single_mut_output(single_mut, search_para, combined_ids, sequence, structure_d):
    
    search_para['pdbid'] = single_mut[0]
    search_para['mut_chain_id'] = single_mut[1]
    search_para['mut_pos'] = single_mut[2]
    
    mut_pos_list = single_mut[2]
    mut_aa_list = single_mut[3]
    DDG = single_mut[4]
    mut_id = single_mut[5]
    
    form = search_para['form']
    
    formalized_pairs = \
        Formalized_contacting(search_para, combined_ids, sequence, structure_d)
    
    # reformat the formalized pairs with more information and transform into wt mut pairs
    Ab_Ag_wt_mut_pairs_list = []
    for one_tuple in formalized_pairs:
        for pos, aa in zip(mut_pos_list, mut_aa_list):
            if one_tuple[0] == pos:
                mut_aa = aa; contact_number = one_tuple[6]
                Ab_Ag_wt_mut_pair = Generate_wt_mut_one(one_tuple, form, mut_aa)
                # Attach the DDG and the mutid
                Ab_Ag_wt_mut_pair.extend([contact_number, DDG, mut_id])
                # Load to the Ab_Ag_wt_mut_pairs_list
                Ab_Ag_wt_mut_pairs_list.append(copy.deepcopy(Ab_Ag_wt_mut_pair))
        
    return Ab_Ag_wt_mut_pairs_list
'''###########################################################################################'''
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
    [[wt_Ab_Ag_pair], [mut_Ab_Ag_pair], contact_number, DDG, mut_id],
    ...]    
'''
 
def Workable_output(workable_input, search_para, combined_ids, sequence, structure_d):
    workable = []

    for one_mutation in workable_input:
        one_mut_workable = []
        for single_mut in one_mutation:
            Ab_Ag_wt_mut_pairs_list =\
                Single_mut_output(single_mut, search_para, combined_ids, sequence, structure_d)
            one_mut_workable.extend(Ab_Ag_wt_mut_pairs_list)
            
        workable.append(copy.deepcopy(one_mut_workable))
    
    return workable
###################################################################################
'''
The input of this module:
    mutation_d: mutation directory, where the mutation information file is located. Refer to 
                the function 'Workable_input' for more details
    sequence: the sequence extracted from the 'AAC_2', The least sequence info it should have
                is the sequence of all the mutation complexes
    combined_ids: a dictionary with keys and values in the form:
                {'pdbid':['H_chains', 'L_chains', 'A_chains']}
    structure_d: the directory where the structure information is located.
    
    search_para: a dictionary, contains the parameters of extrating wild type and mutation pairs.
                Refer to the explaination in related functions
Output:
    workable_output: a list with element list in the following form:
        [..., [wt_Ab_Ag_pair, mut_Ab_Ag_pair, contact_number, DDG, mut_id]]
    each element list contains all the mutaion information related to one measured affinity.
    
    'wt_Ab_Ag_pair' and 'mut_Ab_Ag_pair' can be directly fed into the RBFN model to generate scores 
    about the propensity of them being a core.
'''
if __name__ == '__main__':  

    preliminary_pred = {}; workable_output = {}
      
    mutation_d = '/home/leo/Documents/Database/Data_Code_Publish/Mutations'
#    mutation_d = '/home/leo/Documents/Database/Pipeline_New/Codes/Results'
    workable_input = Workable_input(mutation_d) # Change the file names inside this function.
    
    os.chdir('/home/leo/Documents/Database/Data_Code_Publish/Structures')
    with open('sequence', 'r') as f:
        sequence = json.load(f)
    with open('combined_ids', 'r') as f:
        combined_ids = json.load(f)
        
    search_para = {}
    search_para['moving'] = True
    search_para['step_size'] = 0.25
    search_para['start_dist'] = 3
    search_para['end_dist'] = 8
    search_para['cut_dist'] = 5.5
    
    for form in ['one', 'multiple', 'flanked']:  
        search_para['form'] = form
        for within_range in [True, False]:
            search_para['within_range'] = within_range
            
            print('Working on: '+ form +'    '+ str(within_range))
            preliminary_pred[form+'_WithinRange_'+str(within_range)] = {}
            # Container
            container = []

            structure_d = '/home/leo/Documents/Database/Data_Code_Publish/Structures/imgt'    
            workable = Workable_output(workable_input, search_para, combined_ids, sequence, structure_d)
            workable_output[form+'_WithinRange_'+str(within_range)] = workable

    saving_d = '/home/leo/Documents/Database/Data_Code_Publish/Codes/Results'
    with open('workable_formated', 'w') as f:
        json.dump(workable_output, f)
#        














