import os
import json
import copy
import sys
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
    if mut_chain_id in combined[0]:
        mut_chain_type = 'H'
    elif mut_chain_id in combined[1]:
        mut_chain_type = 'L'
    elif mut_chain_id in combined[2]:
        mut_chain_type = 'A'
        
    # Find the opposite chain_ids and pos
    if mut_chain_type in ['H', 'L']:
        op_chain_ids = combined[2]
    elif mut_chain_type == 'A':
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
            
    return mut_cdn, mut_chain_type, op_cdn
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
            six_c = copy.deepcopy(six)
            six_c = list(six_c)
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
    h_range.append(list(range(50, 71)))
    h_range.append(list(range(99, 129)))
#    h_range = [[25, 37], [50, 71], [99, 129]]
    within_range = search_para['within_range']
    
    if within_range and formalized_pairs != []: 
        for pair in formalized_pairs:               
            if pair[7] == 'H' and pair[0] not in h_range:
                return []
            if pair[7] == 'L' and pair[0] not in l_range:
                return []
            if pair[8] == 'H' and pair[0] not in h_range:
                return []
            if pair[8] == 'L' and pair[0] not in l_range:
                return []

    return formalized_pairs
'''##########################################################################3'''
def Form_one(contact_pairs):
    form_one = []
    if contact_pairs !=[]:
        contact_pairs.sort(key=lambda x : x[6])
        form_one = [contact_pairs[-1]]
    return form_one

def Form_multiple(contact_pairs):
    return contact_pairs

def Form_flanked(contact_pairs, sequence, pdbid):
    flanked = []
    
    form = Form_one(contact_pairs)
    
    if form != []:
        form_one = form[0]
        
        mut_pos = form_one[0]
        mut_chain_id = form_one[2]
        op_pos = form_one[3]
        op_chain_id = form_one[5]
        
        # Find the amino acids before and after
        mut_before_pos = max(0, mut_pos - 1)
        op_before_pos = max(0, op_pos - 1)
        
        flanked_mut = []
        for m in range(mut_before_pos, mut_pos+2):
            flanked_mut.append(sequence[pdbid][mut_chain_id][m])
        
        flanked_op = []
        for o in range(op_before_pos, op_pos + 2):
            flanked_op.append(sequence[pdbid][op_chain_id][o])
            
        # Replace the mut_aa and op_aa in form one with the flanked
        flanked = list(form_one)
        flanked[1] = flanked_mut
        flanked[4] = flanked_op
        flanked = [tuple(flanked)]# Change back into the same form
    
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
            earch_para['mut_pos']
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
    mut_cdn, mut_chain_type, op_cdn = \
        Extract_cdn(pdbid, mut_pos, mut_chain_id, combined_ids, sequence, structure_d)
        
    contact_pairs =\
       Extract_contacting_pairs(mut_cdn, mut_chain_type, op_cdn, search_para)
   
    formalized_pairs =\
       Formalize_contact_pairs(pdbid, sequence, contact_pairs, search_para, combined_ids)
    
    return formalized_pairs
'''#############################################################################'''
# A LITTLE TEST
#search_para = {}
#search_para['moving'] = False
#search_para['step_size'] = 0.25
#search_para['start_dist'] = 5
#search_para['end_dist'] = 8
#search_para['cut_dist'] = 8
#search_para['form'] = 'multiple'
#search_para['within_range'] = False
#search_para['pdbid'] = '1dvf'
#search_para['mut_pos'] = [32]
#search_para['mut_chain_id'] = 'D'
#
#structure_d = '/home/leo/Documents/Database/Data_Code_Publish/Structures/imgt'
#
#os.chdir('/home/leo/Documents/Database/Data_Code_Publish/Structures')
#with open('combined_ids', 'r') as f:
#    combined_ids = json.load(f)
#with open('sequence', 'r') as f:
#    sequence = json.load(f)
#
#formalized_pairs = Formalized_contacting(search_para, combined_ids, sequence, structure_d)




