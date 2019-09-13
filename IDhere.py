import os
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
                    if temp[i] == 'NA':
                        temp[i] =''                       
                if a[0].strip() in id_dict:
                    id_dict[a[0].strip()].append (temp)
                else:
                    id_dict[a[0].strip()] = [temp]                    
    return id_dict 
#########################################################################           
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
#################################################################             
'''
Here_iddict_combineddict, the purpose of this function is to generate the id_dict
       and combined_chain_id_dict for the pdb files in current working directory. 
       This function can be used for checking whether the codes are correct.
Input:
     id_dict, the output of the above function Id_dict
     combined_chain_id_dict, the output of the above function Combined_chain_id_dict
Output:
    here_id_dict, the id_dict for the pdb files in current working directory.
    here_combined_dict, the combined_chain_id_dict for the pdb files in 
    current working directory.
    
'''

def Here_iddict_combineddict(id_dict, combined_chain_id_dict):
    here_id_dict = {}
    here_combined_dict = {}
    names = os.listdir()
    for f in names:
        if len(f) == 8 and f[5:8] == 'pdb':
            if f[:4] in id_dict:
                here_id_dict[f[:4]] = id_dict[f[:4] ]
            if f[:4] in combined_chain_id_dict:
                here_combined_dict[f[:4]] = combined_chain_id_dict[f[:4] ] 
    return here_id_dict, here_combined_dict
##################################################################

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
def witch_hunt(conbined_chain_id, id_dict):
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
# AAC analysis

#len(witches)
#witches
#good_combined_ids
#good_matched_ids
########################################################################
def main(working_directory):
    import os
    os.chdir(working_directory)
    #os.getcwd()
    summary = open('summary.tsv', 'r')# case sensitive
    file = summary.readlines()        
    summary.close
    # take a look at the results
    id_dict = Id_dict (file)
    #len(id_dict)
    #id_dict['4ydl']
    combined_chain_id_dict = Combined_chain_id_dict (id_dict)
    #len(combined_chain_id_dict)
    #combined_chain_id_dict['4ydl']
    here_id_dict, here_conbined_chain_id = Here_iddict_combineddict(id_dict, combined_chain_id_dict)
    #id_dict.keys()
    #here_id_dict.keys()
    witches, good_combined_ids, good_matched_ids = witch_hunt(combined_chain_id_dict, id_dict)
    
    return good_combined_ids, good_matched_ids

#good_combined_ids, good_matched_ids = main('/home/leo/Documents/Database/Pipeline_New/Latest')
#import os
#import json
#os.chdir('/home/leo/Documents/Database/Pipeline_New/Latest')
#with open('good_combined_ids', 'w') as f:
#    json.dump(good_combined_ids, f)
#with open('good_matched_ids', 'w') as f:
#    json.dump(good_matched_ids, f)


#if __name__ == '__main__':
#    training_directory = '/home/leo/Documents/Database/Pipeline/All with peptide 5+ resolution 4A'
#    testing_directory = '/home/leo/Documents/Database/Pipeline/Complexes with Affinity'
#    
#    training_combined_ids, training_matched_ids = main(training_directory)
#    testing_combined_ids, testing_mathed_ids = main(testing_directory)
#    
#    print('The number of all ids is ', len(training_combined_ids), '\n')
#    for key in testing_combined_ids:
#        del training_combined_ids[key]
#        del training_matched_ids[key]
#    print('The number of all training ids is  ', len(training_combined_ids),'\n')
#    print('The number of all testing ids is  ', len(testing_combined_ids), '\n')
#    
#    print('Show the first six training ids\n' )
#    keys = list(training_combined_ids.keys())
#    for key in keys[:6]:
#        print(key)
#        print(training_combined_ids[key])
#        print(training_matched_ids[key])
#    
#    print('\n Show the first six testing ids\n')
#    keys = list(testing_combined_ids.keys())
#    for key in keys[:6]:
#        print(key)
#        print(testing_combined_ids[key])
#        print(testing_mathed_ids[key])
#        
#    # save the results
#    import os
#    import json
#    os.chdir(training_directory)
#    with open('good_combined_ids', 'w') as f:
#        json.dump(training_combined_ids)
#    with open('good_matched_ids', 'w') as f:
#        json.dump(training_matched_ids, f)
#    print('The the training ids is saved to \n', training_directory)
#    
#    os.chdir(testing_directory)
#    with open('good_combined_ids', 'w') as f :
#        json.dump(testing_combined_ids, f)
#    with open('good_matched_ids','w') as f:
#        json.dump(testing_mathed_ids, f)
#    print('The testing ids is saved to \n', testing_directory)
    



