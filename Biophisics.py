'''
This file is to find the structures of the complexes which can be used to 
explain the biochemical properties of the cores
'''
import os, json, copy
from pyexcel_ods import get_data, save_data

os.chdir('/home/leo/Documents/Database/Pipeline_New/Complexes/Cores')
with open('training_1_1_0_0_1_2_1perchain', 'r') as f:
    training_1_1 = json.load(f)
    
os.chdir('/home/leo/Documents/Database/Pipeline_New/Complexes/Results')
with open('one_to_one_frequency', 'r') as f:
    one_to_one_frequency = json.load(f)

one_to_one_frequency[:3]
training_1_1

'''
Find:
    this function is to find the original core information of the top_n 
    cores, which could be used to find the biophysical informations in the structures
    of the complexes.
Input:
    training_1_1:
        the cores extracted with match_type 1_1
    one_to_one_frequency:
        a list, gives the frequencies of the cores with the match_type 1_1
    top_n:
        an integer, gives how many of the top frequency cores will be considered.
Output:
    core_partition:
        a dictionary, with the keys '1', '2', ..., which give the rank of the frequencies
         and the values the lists of the extracted cores
    
'''
def Find(training_1_1, one_to_one_frequency, top_n = 3):
    core_partition ={}
    for i, fre_core in enumerate(one_to_one_frequency[:top_n]):
        temp_container = []
        for core in training_1_1:
            if core[0][0] == fre_core[0] and core[1][0] == fre_core[1]:
                temp_container.append(core)
        temp_container.sort(key = lambda x:x[2], reverse=True)
        core_partition[str(i+1)]=copy.deepcopy(temp_container)
            
    return core_partition

#core_partition = Find(training_1_1, one_to_one_frequency, top_n = 3)
#core_partition.keys()
#len(core_partition['1'])
#len(core_partition['2'])
#core_partition['1']
#os.chdir('/home/leo/Documents/Database/Pipeline_New/Complexes')
#with open('contact', 'r') as f:
#    contact = json.load(f)

#############################################################################
'''
Find the mutations with the number of mutated amino acids no less than 2
'''
def Find_more_than_1_continuous_mutation():
    os.chdir('/home/leo/Documents/Database/Pipeline_New/Codes/Results')
    data_raw = get_data('Keating.ods')
    keys = data_raw.keys()
    data_multimutation = {}  
    end = 'affinities'  
    for key in keys:
        data_multimutation[key] = []
        for ind, mut in enumerate(data_raw[key]):
            if mut[0] == key:
                start_n = ind
            elif mut[0] == end:
                end_n = ind
                if end_n - start_n >= 2:
                    target_set = data_raw[key][start_n:end_n]
                    # separate the target_set into different mutations on different chains
                    chain_names = []
                    for ts in target_set:
                        if ts[1] not in chain_names:
                            chain_names.append(ts[1])
                    # set a count dictionary and count the number of mutations in each chain
                    count ={}
                    for name in chain_names:
                        count[name] = []
                    for ts in target_set:
                        if ts[1] in count:
                            count[ts[1]].append(ts[2])
                    # Find the smallest length in the count dictionary
                    continuous = {}
                    for ke, value in count.items():
                        continuous[ke] = 0
                        value.sort()
                        start_index = value[0]
                        for index in value:
                            if index - start_index == 1:
                                continuous[ke] = 1
    
                    # Add those with all the values in the continuous be 1
                    add = True
                    for c, v in continuous.items():
                        if v == 0:
                            add = False
                    if add:
                        data_multimutation[key].extend(data_raw[key][start_n:end_n+1])
    delete = []                    
    for pdbid, mut in data_multimutation.items():
        if mut == []:
            delete.append(pdbid)
    for d in delete:
        del data_multimutation[d]
    
    return data_multimutation        

data_multimutation = Find_more_than_1_continuous_mutation()
data_multimutation                
# Save the results as ods
#os.chdir('/home/leo/Documents/Database/Pipeline_New/Codes/Results')
#save_data('multimutations.ods', data_multimutation)
          

  
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    