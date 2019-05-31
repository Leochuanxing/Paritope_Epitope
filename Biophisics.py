'''
This file is to find the structures of the complexes which can be used to 
explain the biochemical properties of the cores
'''
import os, json, copy

os.chdir('/home/leo/Documents/Database/Pipeline_New/Complexes/Cores')
with open('training_1_1_0_0_1_2_1perchain', 'r') as f:
    training_1_1 = json.load(f)
    
os.chdir('/home/leo/Documents/Database/Pipeline_New/Complexes/Results')
with open('one_to_one_frequency', 'r') as f:
    one_to_one_frequency = json.load(f)

one_to_one_frequency
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

core_partition = Find(training_1_1, one_to_one_frequency, top_n = 3)
core_partition
len(core_partition['1'])
