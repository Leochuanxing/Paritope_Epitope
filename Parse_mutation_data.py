'''
REMARK:
    1. This file is to parse the mutation data given by the paper AB-Bind: Antibody binding mutational
        database for computational affinity predictions
    2. The complex with 1dvf is a complex between two antibodies. The parsing of the mutations on this
        complex was done manually.
'''
import json
import xlrd
import os
import copy
from pyexcel_ods import get_data, save_data
os.chdir('/home/leo/Documents/Database/Data_Code_Publish/Codes')
from AAC_2 import Chain_seq
#######################################################################
def T_S(triple):
    TripleSingle =  [['TYR', 'Y'], ['LYS', 'K'],['ASP', 'D'], ['ASN', 'N'], ['TRP', 'W'], ['PHE', 'F'], ['GLN', 'Q'],
       ['GLU', 'E'], ['PRO', 'P'], ['GLY', 'G'], ['THR', 'T'],['SER', 'S'], ['ARG', 'R'], ['HIS', 'H'],
       ['LEU', 'L'], ['ILE', 'I'], ['CYS', 'C'], ['ALA', 'A'], ['MET', 'M'], ['VAL', 'V']]
    for TS in TripleSingle:
        if TS[0] == triple:
            return TS[1]
def S_T(single):
    TripleSingle =  [['TYR', 'Y'], ['LYS', 'K'],['ASP', 'D'], ['ASN', 'N'], ['TRP', 'W'], ['PHE', 'F'], ['GLN', 'Q'],
       ['GLU', 'E'], ['PRO', 'P'], ['GLY', 'G'], ['THR', 'T'],['SER', 'S'], ['ARG', 'R'], ['HIS', 'H'],
       ['LEU', 'L'], ['ILE', 'I'], ['CYS', 'C'], ['ALA', 'A'], ['MET', 'M'], ['VAL', 'V']]
    for TS in TripleSingle:
        if TS[1] == single:
            return TS[0]
##################################################################################

#############################################################################33
'''
The mutations positions are give as the original pdb files, for example, some of the 
chains may not begin with 1. Therefore we have to index the positions  of the mutaions with 
the positions in our analysis

Match_index:
    a function to match the indices in the mutation data with the mutations in our system
Input:
    file: a pdbfile
Output:
    match_indices_dict:
    a dictionary, with keys the chain names and the values in the following form
    [[101, '100a'], [102, '100b'], ...]
'''

def Match_index(file, combined_chain_id):
    match_indices_dict = {}
    # creat a tracker dictionary, and a counter dictionary
    tracker = {}
    counter = {}
    ids = combined_chain_id[0] + combined_chain_id[1] + combined_chain_id[2]
    for i in ids:
        tracker[i] = ''
        counter[i] = -1
        match_indices_dict[i] = []
    for line in file:
        if line[:6] == 'ATOM  ' and line[21] in ids:
            # update the parameters
            if tracker[line[21]]!= line[22:27].strip():
                counter[line[21]] += 1
                tracker[line[21]] = line[22:27].strip()
                count = counter[line[21]] 
                track = tracker[line[21]]
            # extract all the parameters corresponding to line[21]
                match_indices_dict [line[21]].append(copy.deepcopy([count, track]))
                                    
    return match_indices_dict
    
#######################################################################
'''
Parse the mutation data so that they can be easily used
Parse:
     A function to parse the mutation data
Input:
    sheet: the mutation data downloaded from the paper mentioned at the beginning, with the 
    non antibody-antigen complex excluded, they are ['1ak4', '1ffw', '1jtg', '1ktz', '1t83', '3k2m', '3wjj']
    and the ones begin with hm
Output:
    parse_dictionary_list:
    a dictionnary with keys pdbids and values a list of dictionaries, of which the keys are
    'mutations' and 'affinities'
    for example
    the value of 'mutations' is in the form [pdbid, chain_id, pos, mutation_aa, wt_aa]

'''

def Parse(sheet, matched_ids):
    parse_dictionary_list=[]
    # Lets collect the mutation ids first
    for i in range(1,1015):
        pdbid = sheet.cell_value(i,0).lower()        
        mut_set = (sheet.cell_value(i,4)).split(',')
        dG = sheet.cell_value(i,5)
        if pdbid in matched_ids:
            temp_dict = {} 
            temp_dict['pdbid'] = pdbid
            temp_dict['mutations'] = []
            temp_dict['mutation_id'] = i
            for mut in mut_set:                
                chain_aa_pos_aa = mut.split(':')
                chain_id = chain_aa_pos_aa[0]        
                wt_aa = S_T(chain_aa_pos_aa[1][0])
                pos = chain_aa_pos_aa[1][1:][:-1]
                mutation_aa = S_T(chain_aa_pos_aa[1][-1])
                temp_dict['mutations'].append([pdbid, chain_id, pos, mutation_aa, wt_aa])
            temp_dict['affinities'] = dG
            parse_dictionary_list.append(copy.deepcopy(temp_dict))
    
    return parse_dictionary_list

###############################################################################
'''
Check_indices:
    A function to check whether the matched indices are correctly matched
Input:
    sequence, pdbid, matched_indices_dict
    *************************
    the matched_indices_dict should be corresponding to the pdbid
Output:
    unmatch:
        a list shows the unmatched indices, if all the indices are matched up, the unmatch
        is an empty list.
'''
def Check_indices(sequence, matched_indices, onepdb_parsed_dictionary):
    unmatch = []
    for single in onepdb_parsed_dictionary:        
        pdbid = single['pdbid']
        chain_sequence = sequence[pdbid]
        mutations = single['mutations']
        # Check
        for mu in mutations:
            data_index = mu[2]
            chain = mu[1]
            for i in matched_indices[chain]:
                if i[1].lower() == data_index:
                    my_index = i[0]
                    
            if mu[4] != chain_sequence[chain][my_index]:
                unmatch.append(mu)
            
    return unmatch
#############################################################################
'''
Direct match:
    a function to use the indices directly from the mutation data without the 
    converting by the matched_indices
'''
def Direct_match(sequence, onepdb_parsed_dictionary):
    unmatch = []
    for single in onepdb_parsed_dictionary:
        pdbid = single['pdbid']
        chain_sequence = sequence[pdbid]
        mutations = single['mutations']
        # Check
        for mu in mutations:
            data_index = int(mu[2])
            chain = mu[1]
            data_aa = mu[4]
            if chain_sequence[chain][data_index] != data_aa:
                unmatch.append(mu)               
    return unmatch

##############################################################################
def Match_Dic(parse_dictionary_list):   
    # Get all the pdbids
    dic_all_mutations = {}
    
    for dic in parse_dictionary_list:
        if dic['pdbid'] not in dic_all_mutations:
            dic_all_mutations[dic['pdbid']] = [dic]
        elif dic['pdbid'] in dic_all_mutations:
            dic_all_mutations[dic['pdbid']].append(dic)
            
    os.chdir('/home/leo/Documents/Database/Data_Code_Publish/Structures/imgt')
    match_indices_all = {}
    for key in dic_all_mutations:        
        with open(key+'.pdb', 'r') as file:    
            combined_chain_id = combined_ids[key]    
            match_indices_all[key] = Match_index(file, combined_chain_id)
            
    return match_indices_all, dic_all_mutations
############################################################################
def Use_match(matched_indices, onepdb_parsed_dictionary):

    for single in onepdb_parsed_dictionary:        
        mutations = single['mutations']
        # Check
        for mu in mutations:
#            data_index = mu[2]
            chain = mu[1]
            for i in matched_indices[chain]:
#                print(mu[2])
#                if mu[2] == '100a':
#                    print(mu[0], i[1], mu[2])
                if i[1].lower() == mu[2]:
#                    print(mu[2])
#                    my_index = i[0]
                    mu[2] = i[0]
    return onepdb_parsed_dictionary
###################################################################
'''move up y n'''
def Use_unmatch(n, onepdb_parsed_dictionary):
    for single in onepdb_parsed_dictionary:
        mutations =  single['mutations']
        # Check different chains
        for mu in mutations:
            mu[2] = int(mu[2])+n
#####################################################################
def Format(onepdb_parsed_dictionary):
    formated = []
    for single in onepdb_parsed_dictionary:
        pdbid = single['pdbid']
        mutations = single['mutations']
        affinities = single['affinities']
        mutation_id = single['mutation_id']
        for i in range(len(mutations)):
            if i ==0:
                mu = mutations[i]
                mu[0] = pdbid
                mu.append(mutation_id)
                formated.append(mu)
            else:
                mu = mutations[i]
                mu[0] = ''
                formated.append(mu)
        formated.append(['affinities', affinities, 1, 'DDG'])
    return formated

#####################################################################
'''Run the following step by step can parse the mutation information'''
'''
Set up a index so that the mutaion data can be referred to the data in my research

Do the following step by step
'''  
          
#unmatch = Check_indices(sequence,  matched_indices, onepdb_parsed_dictionary)
#unmatch
#sequence[pdb]['I'][94]
#len(sequence[pdb]['I'])
def Change_ind_to_match(sequence, matched_indices, onepdb_parsed):
    onepdb_parsed_dictionary = copy.deepcopy(onepdb_parsed)
    onepdb_parsed_dictionary = Use_match(matched_indices, onepdb_parsed_dictionary)
    #onepdb_parsed_dictionary = Use_unmatch(onepdb_parsed_dictionary)
#    onepdb_parsed_dictionary
    
    direct_unmatch = Direct_match(sequence, onepdb_parsed_dictionary)

    n = -1
    while direct_unmatch !=[] and n <= 1:
        
        onepdb_parsed_dictionary = copy.deepcopy(onepdb_parsed)
        onepdb_parsed_dictionary = Use_match(matched_indices, onepdb_parsed_dictionary)
        Use_unmatch(n, onepdb_parsed_dictionary)
        direct_unmatch = Direct_match(sequence, onepdb_parsed_dictionary)
#        onepdb_parsed_dictionary = Use_match(matched_indices, onepdb_parsed_dictionary)
        n += 1
        
    return direct_unmatch, onepdb_parsed_dictionary

def Change_chain_name(sequence, combined_id_one, matched_indices, onepdb_parsed):
    Hchain = combined_id_one[0]
    Lchain = combined_id_one[1]
    Achain = combined_id_one[2]
    for single in onepdb_parsed:
        mutations = single['mutations']
        for mu in mutations:
            if mu[1] == 'H':
                mu[1] = Hchain[0]
            elif mu[1] in 'L':
                mu[1] = Lchain[0]  
            elif mu[1] in 'A':
                mu[1] = Achain[0]
    direct_unmatch, onepdb_parsed_dictionary = Change_ind_to_match(sequence, matched_indices, onepdb_parsed)
    
    return direct_unmatch, onepdb_parsed_dictionary
###########################################################################################
def Format_affinities(data):
    for key, value in data.items():
        for mu in value:
            if mu != [] and mu[0] == 'affinities' and mu[3] == 'DDG':
                mu[2] = ''
#                mu[1] = 0
#                mu[2] = dg
'''##################################################################################'''
if __name__ == '__main__':
    
    loc = '/home/leo/Documents/Database/Data_Code_Publish/Mutations/Mutation.xlsx'
    wb = xlrd.open_workbook(loc) 
    sheet = wb.sheet_by_index(0) 
    # Update the data with 1dvf  
    os.chdir('/home/leo/Documents/Database/Data_Code_Publish/Structures')
    with open('matched_ids', 'r') as f:
        matched_ids = json.load(f)
    with open('combined_ids', 'r') as f:
        combined_ids = json.load(f)
    with open('sequence', 'r') as f:
        sequence = json.load(f)
        
    matched_1dvf = {}
    matched_1dvf['1dvf'] = [['B', 'A', 'C'], ['B', 'A', 'D']]
    combined_1dvf = {}
    combined_1dvf['1dvf'] = ['B', 'A', 'CD']
    # Get the sequence of idvf
    os.chdir('/home/leo/Documents/Database/Data_Code_Publish/Structures/imgt')
    with open('1dvf.pdb', 'r') as file:
        sequence_1dvf = Chain_seq(file, combined_1dvf['1dvf'])
    sequence_1dvf_dict = {}
    sequence_1dvf_dict['1dvf'] = sequence_1dvf
    
    matched_ids.update(matched_1dvf)
    combined_ids.update(combined_1dvf)
    sequence.update(sequence_1dvf_dict)

    # Check if the sequences matched 
    os.chdir('/home/leo/Documents/Database/Data_Code_Publish/Mutations') 
    parse_dictionary_list = Parse(sheet, matched_ids)
#    parse_dictionary_list[:6]   
    match_indices_dict, dic_all_mutations = Match_Dic(parse_dictionary_list)
#    len(match_indices_dict)
#    len(dic_all_mutations)
    keys = list(match_indices_dict.keys())

    # Find the final qualified mutaions
    pdb_qualified = []; parsed_dictionaries = []
    for pdb in keys:
#        pdb = keys[6]                
        onepdb_parsed = dic_all_mutations[pdb]
        matched_indices = match_indices_dict[pdb] 
        combined_id_one = combined_ids[pdb]
        
        try:
            direct_unmatch, onepdb_parsed_dictionary = Change_ind_to_match(sequence, matched_indices, onepdb_parsed)
        except KeyError:
            try:
                direct_unmatch, onepdb_parsed_dictionary = Change_chain_name(sequence, combined_id_one, matched_indices, onepdb_parsed)
            except :
                direct_unmatch = ['Not good']
                print(pdb)
        except :
            direct_unmatch = ['Not good']
            print(pdb)
        
        if direct_unmatch == []:
            pdb_qualified.append(pdb)
            parsed_dictionaries.append(onepdb_parsed_dictionary)
        else:
            print(pdb)
            
    # Format the dictionary before saving
    formated_dictionary = {}
    for one_pdb in parsed_dictionaries:
        formated = Format(one_pdb)
        formated_dictionary[one_pdb[0]['pdbid']] = copy.deepcopy(formated)
    # Change the format of the affinity a little bit   
    Format_affinities(formated_dictionary)
    # Save the results    
    os.chdir('/home/leo/Documents/Database/Data_Code_Publish/Mutations')   
    save_data('Formated.ods', formated_dictionary)
    
#formated_dictionary
#########################################################################

    
