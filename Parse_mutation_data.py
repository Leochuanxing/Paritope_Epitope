'''
This file is to parse the mutation data given by the paper AB-Bind: Antibody binding mutational
database for computational affinity predictions
'''
import json
import xlrd
import os
import copy
from pyexcel_ods import get_data, save_data
loc = '/home/leo/Documents/Database/Pipeline_New/Mutations/Mutation data.xlsx'
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

#os.chdir('/home/leo/Documents/Database/Pipeline_New/Mutations')
#data_raw = get_data('Mutations.ods')
#data = []
#for d in data_raw:
#    if d != []:
#        data.append(d)
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
loc = '/home/leo/Documents/Database/Pipeline_New/Mutations/Mutation.xlsx'
wb = xlrd.open_workbook(loc) 
sheet = wb.sheet_by_index(0) 

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
Check if the indices are correct

Check_indices:
    A function to check whether the indices are correct
Input:
    sequence, pdbid, matched_indices_dict
    *************************
    the matched_indices_dict should be corresponding to the pdbid
Output:
    unmatch:
        a list shows the unmatched indices, if all the indices are matched up, the unmatch
        is an empty list.
'''
def Check_indices(sequence,  matched_indices, onepdb_parsed_dictionary):
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

 

#########################################################################
os.chdir('/home/leo/Documents/Database/Pipeline_New/All with peptide 5+ resolution 4A')
with open('sequence', 'r') as f:
    sequence = json.load(f)
with open('combined_ids', 'r') as f:
    combined_ids = json.load(f)
##############################################################################
def Match_Dic(parse_dictionary_list):   
    # Get all the pdbids
    dic_all_mutations = {}
    
    for dic in parse_dictionary_list:
        if dic['pdbid'] not in dic_all_mutations:
            dic_all_mutations[dic['pdbid']] = [dic]
        elif dic['pdbid'] in dic_all_mutations:
            dic_all_mutations[dic['pdbid']].append(dic)
            
    os.chdir('/home/leo/Documents/Database/Pipeline_New/All with peptide 5+ resolution 4A/structure')
    match_indices_all = {}
    for key in dic_all_mutations:        
        with open(key+'.pdb', 'r') as file:    
            combined_chain_id = combined_ids[key]    
            match_indices_all[key] = Match_index(file, combined_chain_id)
            
    return match_indices_all, dic_all_mutations
############################################################################
#def My_pdbid():
#    total_number_mutation = 0
#    pdbids = []
#    os.chdir('/home/leo/Documents/Database/Pipeline_New/Codes/Results')
#    data = []
#    for i in [1,2,3]:
#        sheet = 'Sheet'+str(i)
#        data_raw = get_data('Affinities.ods')[sheet]    
#        for d in data_raw:
#            if d != []:
#                data.append(d)
#    for mu in data:
#        if mu[0] != '' and mu[0] != 'affinities':
#            total_number_mutation += 1
#            if mu[0] not in pdbids:
#                pdbids.append(mu[0])
#    
#    return total_number_mutation, pdbids
######################################################################
#def Pdbid_in_other_data(my_ids, other_ids):
#    more_ids = []
#    for i in other_ids:
#        if i not in my_ids:
#            more_ids.append(i)
#    return more_ids
########################################################################
def Use_match(matched_indices, onepdb_parsed_dictionary):

    for single in onepdb_parsed_dictionary:        
        mutations = single['mutations']
        # Check
        for mu in mutations:
            data_index = mu[2]
            chain = mu[1]
            for i in matched_indices[chain]:
                if i[1].lower() == data_index:
                    my_index = i[0]
                    mu[2] = my_index
    return onepdb_parsed_dictionary
###################################################################
def Use_unmatch(onepdb_parsed_dictionary):
    for single in onepdb_parsed_dictionary:
        mutations =  single['mutations']
        for mu in mutations:
            mu[2] = int(mu[2])-1
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
#    os.chdir('/home/leo/Documents/Database/Pipeline_New/Results')   
#    data = get_data('Affinities.ods')['Sheet2']
#    data.update({"Sheet 4": formated})
#    save_data('Affinities.ods', data)
#####################################################################
def Change_chain_name(onepdb_parsed_dictionary):
    for single in onepdb_parsed_dictionary:
        mutations = single['mutations']
        for mu in mutations:
            if mu[1] == 'A':
                mu[1] = 'C'
            if mu[1] == 'H':
                mu[1] = 'B'
        
###################################################################
'''Run the following step by step can parse the mutation information'''
'''
Set up a index so that the mutaion data can be referred to the data in my research

Do the following step by step
'''  
os.chdir('/home/leo/Documents/Database/Pipeline_New/All with peptide 5+ resolution 4A')
with open('matched_ids', 'r') as f:
    matched_ids = json.load(f)
os.chdir('/home/leo/Documents/Database/Pipeline_New/Mutations') 
parse_dictionary_list = Parse(sheet, matched_ids)
parse_dictionary_list[:6]             
 
match_indices_dict, dic_all_mutations = Match_Dic(parse_dictionary_list)
len(match_indices_dict)
len(dic_all_mutations)
keys = list(match_indices_dict.keys())
keys
#
#
pdb = keys[3]
pdb
onepdb_parsed_dictionary = dic_all_mutations[pdb]
matched_indices= match_indices_dict[pdb] 
onepdb_parsed_dictionary
#matched_indices
#Change_chain_name(onepdb_parsed_dictionary)
combined_ids[pdb]
unmatch = Check_indices(sequence,  matched_indices, onepdb_parsed_dictionary)
unmatch
onepdb_parsed_dictionary = Use_match(matched_indices, onepdb_parsed_dictionary)
#onepdb_parsed_dictionary = Use_unmatch(onepdb_parsed_dictionary)
#
direct_unmatch = Direct_match(sequence, onepdb_parsed_dictionary)
direct_unmatch
Use_unmatch(onepdb_parsed_dictionary)
#
#
formated = Format(onepdb_parsed_dictionary)
##formated
#
#os.chdir('/home/leo/Documents/Database/Pipeline_New/Codes/Results')   
#data = get_data('Keating.ods')
#len(data)
#data.update({pdb: formated})
#save_data('Keating.ods', data)


#onepdb_parsed_dictionary
#sequence[pdb]['B'][97]
#matched_indices['I']
#matched_indices.keys()


#########################################################################
'''*****The following code is to Format the affinities*********'''
#os.chdir('/home/leo/Documents/Database/Pipeline_New/Codes/Results')
#data = get_data('Keating.ods')
#len(data)
#data.keys()
def Format_affinities(data):
    for key, value in data.items():
        for mu in value:
            if mu != [] and mu[0] == 'affinities' and mu[3] == 'DDG':
                dg = mu[1]
                mu[1] = 0
                mu[2] = dg
#Format_affinities(data)   
#data
#save_data('Keating.ods', data)      
