'''
This file is to extract the sequences of the chains in the complex, the coordinates of 
the atoms in the antigen chain and the coordinates of the atoms of the CDRs.
'''
###################################################
import json    
import os
import math
#########################################################
#os.chdir('/home/leo/Documents/Database/Pipeline/All with peptide 5+ resolution 4A')
#file = open('summary.tsv', 'r')
#summary = file.readlines()
#file.close
################################################
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
       
def Id_dict (summary):
    id_dict = {}
    for line in summary:        
        vector = line.split('\t')
#        Make sure the vector is long enough
        if len(vector) >= 10:            
#        Deal with the | in a[4]            
            for i in vector[4].split('|') :
                temp = [vector[1].strip(), vector[2].strip(), i.strip()]
                for i in range(0,3):
                    if temp[i] == 'NA':
                        temp[i] =''                       
                if vector[0].strip() in id_dict:
                    id_dict[vector[0].strip()].append (temp)
                else:
                    id_dict[vector[0].strip()] = [temp]   
    # Remove the item with key 'pdb'
    id_dict.pop('pdb')
    return id_dict  
##############################################
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
#####################################################################
'''
witch_hunt,. is to find out the pathological complex with the same chain shows up
          both as the antigen chain and the antibody chain.
inputs: conbined_chain_id, a dictionary with keys pdb id, and values in the form
        of ['BDF', 'ACE', 'GH']
return: witches, a list of weird pdb ids
        GoodPeople, a dictionary with keys of good pdb ids, and values in the form
        of ['BDF', 'ACE', 'GH']
'''
def witch_hunt(combined_chain_id, id_dict):
    # creat a container to contain the witches and GoodPeople
    # the witches are the ones with the same chain appears both as the antigen chain
    # and as the antibody chain.
    witches = []
    for pdbid in combined_chain_id:
        antibody_ids = combined_chain_id[pdbid][0] + combined_chain_id[pdbid][1] 
        antigen_ids = combined_chain_id[pdbid][2]
        witch = False
        # here is to find the witch
        for chain in antibody_ids:
            if chain in antigen_ids:
                witch = True
                break
        if witch:
            witches.append(pdbid)
    for witch in witches:
        del combined_chain_id[witch]
        del id_dict[witch]
        
            
    return combined_chain_id, id_dict
###########################################################################
'''
Keep_the_latest:
    a function to select the latest complexes 
Input:
    summary:
        a summary file corresponding to the complexes
    percent:
        gives the percentage of the latest complexes
Output:
    training_ids:
        a list contains the training ids
    testing_ids:
        a list contains the testing ids
'''
def Keep_the_latest(summary, percent):
    pdbid_date_dict = {}
    for line in summary:        
        vector = line.split('\t')
#        Make sure the vector is long enough
        if len(vector) >= 10 and vector[0] != 'pdb': 
            # Change the date into number form 
            date_str = vector[9]
            vector_date_str = date_str.split('/')
            year = float(vector_date_str[2])
            if year <= 30:
                year += 2000
            else:
                year += 1900
                
            month = float(vector_date_str[0])
            day = float(vector_date_str[1])
            
            date_number = 10000*year+100*month+day
            
            if vector[0] not in pdbid_date_dict:
                pdbid_date_dict[vector[0]] = date_number
            
    # select the training ids and the testing ids from the pdbid_date_dict
    pdbid_date_list = []
    for key, value in pdbid_date_dict.items():
        pdbid_date_list.append([key, value])
    pdbid_date_list.sort(key = lambda x:x[1])
    
    # Calulate the number of training set
    number_training_ids = math.floor(len(pdbid_date_list)*(1-percent))
    # load up the training_ids and the testing_ids
    training_ids = []
    testing_ids = []
    for i in range(len(pdbid_date_list)):
        if i <= number_training_ids:
            training_ids.append(pdbid_date_list[i][0])
        else:
            testing_ids.append(pdbid_date_list[i][0])
            
    return training_ids, testing_ids
################################################
# Extract the chain sequences 
'''
Chain_seq is to extract sequences for all the chains in the complex
Inputs 
        file, a pdb file
        combined_chain_id, a list of the form ['BDF', 'ACE', 'GH']
         in the order of heavy chains, light chains, and antigen chains
Returns: seq, a dictionary of sequences, with the chain id as keys
'''
def Chain_seq(file, combined_chain_id):
    # Combine all the ids together
    ids = combined_chain_id[0] + combined_chain_id[1] + combined_chain_id[2]
    # creat an empty dictionary, set a tracker to track whether an aa should be
    # added to the seq of a partitular chain
    seq = {}
    tracker = {}
    for i in ids:
        seq[i] = []
        tracker[i] = ''
    # load the sequences
    for line in file:
        if line[:6] == 'ATOM  ' and line[21] in ids:
            """Make sure only record the aa when the position number changes"""
            if tracker[line[21]] != line[22:27]:
                seq[line[21]].append(line[17:20])
                tracker[line[21]] = line[22:27]
    return seq


#with open('sequence', 'w') as f:
#    json.dump(sequence, f)

#with open('sequence', 'r') as f:
#    sequence = json.load(f)
#len(sequence)
#sequence['4cad'].keys()
#len(sequence['1adq']['A'])
###############################################################

'''
Coordinates is to extract the coordinates used to calculate the interactions.
inputs: file, a pdb file
        id_dict, combined_chain_id, a list of the form ['BDF', 'ACE', 'GH']
        in the order of heavy chains, light chains, and antigen chains
return: cdn, a dictionary in the form of with keys ['h1H', 'h2H', 'h3H',
         'l1L', 'l1L', 'l1L', ..Antigen chain ids..]
        and the coordinates are in the form of [15.1, 2.2, 3.2, pos, aa]
        pos is an integer, indicates the position in the corresponding chain.
        aa, is the name of the amino acid.
'''
def Coordinates(file, combined_chain_id):
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
    
    # set the range of CDRh and CDRl, all the numbers take the same counting system
    # as python, with the firt one numbered 0. The following method is a conserved way.
    # The starting and the ending point covers all the cdrs nomatter in what way we 
    # define the cdrs.
    #fluctuate the lower limit and upper limit by 3
    l_range = [[23, 40], [49, 63], [89, 110]]
    h_range = [[25, 37], [50, 71], [99, 129]]
#    l_range = [[20, 43], [46, 66], [86, 110]]
#    h_range = [[22, 40], [47, 74], [96, 132]]
    
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
                #Tell the CDR type and load
                if count in range(h_range[0][0], h_range[0][1]+1):
                    cdn['h1'+line[21]].append([float(line[30:38]), float(line[38:46]), float(line[46:54]),
                                       count, line[17:20]])
                if count in range(h_range[1][0], h_range[1][1]+1):
                    cdn['h2'+line[21]].append([float(line[30:38]), float(line[38:46]), float(line[46:54]),
                                       count, line[17:20]])
                if count in range(h_range[2][0], h_range[2][1]+1):
                    cdn['h3'+line[21]].append([float(line[30:38]), float(line[38:46]), float(line[46:54]),
                                       count, line[17:20]])
    
            if c_type == 'L':
                #Tell the CDR type and load
                if count in range(l_range[0][0], l_range[0][1]+1):
                    cdn['l1'+line[21]].append([float(line[30:38]), float(line[38:46]), float(line[46:54]),
                                       count, line[17:20]])
                if count in range(l_range[1][0], l_range[1][1]+1):
                    cdn['l2'+line[21]].append([float(line[30:38]), float(line[38:46]), float(line[46:54]),
                                       count, line[17:20]])
                if count in range(l_range[2][0], l_range[2][1]+1):
                    cdn['l3'+line[21]].append([float(line[30:38]), float(line[38:46]), float(line[46:54]),
                                       count, line[17:20]])
            if c_type == 'A':
                cdn[line[21]].append([float(line[30:38]), float(line[38:46]), float(line[46:54]),
                                       count, line[17:20]])            
    
    return cdn
'''Extract all the coordinates and store them in dictionary coordinates with keys
   pdbid and the elements cdn
'''
#good_combined_ids['2hh0']

'''Take a look at the results'''

#coordinates['1adq'].keys()
#coordinates['1adq']['h1H']
##############################################################
# Extract the contact
'''
Get_contact is to calculate the contact number between amino acids. This function
        is time consuming.
inputs, cdn, a dictionary in the form of with keys ['h1H', 'h2H', 'h3H',
         'l1L', 'l1L', 'l1L', ..Antigen chain ids..] 
         and the coordinates are in the form of [15.1, 2.2, 3.2, pos, aa]
        pos is an integer, indicates the position in the corresponding chain
        aa, is the name of the amino acid.
        cutoff, a float, gives the cutoff distance
        matched_ids, a list in the form of [[H,L,A], [L, M, N]], where 
        [H, L, A] means those three are in a contacting group
return: contact, a list, in the form of [[h1HA, 32, 15, 8], ....]
        this means, the amino acid at position 32, which is located at CDRh1 of chain H, 
        contact with amino acid at position 15 of chain 'A'. The contact number 
        under the given cutoff is 8. The contact number is calculated by the following way:
            if atomA1 from aaA contacts with atomB1 from aaB, then the contact number 
            between aaA and aaB increased by 1. The contact between atomA1 and atomB1
            is only counted once.
'''
def Get_contact(cdn, matched_ids, cutoff = 4):
    # Creat an empty list to contain the temporary results
    contact_temp = []
    squared_cutoff = cutoff**2
    # sorting the keys into CDR and antigen groups
    # it is better to use the information of the matched ids
    # the grouped results should be stored in the form of[ [[h1H, h2H,h3H], [A]], ...]
    grouped =[]
    for matched in matched_ids:
        if matched[2] != '':
            if matched[0] != '':
                grouped.append([['h1'+matched[0], 'h2'+matched[0], 'h3'+matched[0]], [matched[2]]])
            if matched[1] != '':
                grouped.append([['l1'+matched[1], 'l2'+matched[1], 'l3'+matched[1]], [matched[2]]])
    #calculate the contact according to the grouped
    for match in grouped: 
    # calculate the distance and iterating through all possible combinations
        for i in match[0]:
            for j in match[1]:
                for atom1 in cdn[i]:
                    for atom2 in cdn[j]:
                        # We can accelerate this process by selecting the max abs first
                        diff = [atom1[0]-atom2[0],atom1[1]-atom2[1],atom1[2]-atom2[2]]                        
                        # is it faster to compare the square than the sequare root?
                        s = 0
                        for component in diff:
                            s += component**2
                            if s > squared_cutoff:# this step can accelerate the calculation by a lot.
                                break                        
                        if s <= squared_cutoff:
                            contact_temp.append([i+j, atom1[3], atom2[3]])         
    # Count method: Creat a dictionary to count\
    contact = []
    count_dict = {}
    for i in contact_temp:
        string = i[0] + '_' + str(i[1]) + '_' + str(i[2])
        if string in count_dict:
            count_dict[string] += 1
        else:
            count_dict[string] = 1
    # change the count_dict to contact
    contact = []
    for i in count_dict:
        element = i.split('_')
        element[1] = int(element[1])
        element[2] = int(element[2])
        element.append(count_dict[i])
        contact.append(element)
            
    return contact
##################################################################

###############################################
'''
wd: the working directory
percent: the percentage of the latest ids
'''
def main(wd, percent):
    # Extract the information from the summary file
    os.chdir(wd)
    file = open('summary.tsv', 'r')
    summary = file.readlines()
    file.close
    
    id_dict = Id_dict (summary)
    combined_chain_id_dict = Combined_chain_id_dict (id_dict)
    combined_ids, matched_ids = witch_hunt(combined_chain_id_dict, id_dict)
    training_ids, testing_ids = Keep_the_latest(summary, percent)    

    # Extract the sequence
    sequence = {}
    for i in combined_ids:
        with open(i + '.pdb', 'r') as file:
            sequence[i] = Chain_seq(file, combined_ids[i]) 
    
    # Extract the coordinates
    coordinates = {}
    for i in combined_ids:
        with open(i + '.pdb', 'r') as file:
            coordinates[i] = Coordinates(file, combined_ids[i])
    
    #Fine the contact between the antigens and the CDRs
    import time
    start =time.clock()
    contact = {}
    n = 0
    for i in matched_ids:
        n += 1
        print('Calculating   ' + i + '     ' + str(n))
        contact[i] = Get_contact(coordinates[i], matched_ids[i], cutoff = 4)
    end = time.clock()
    print('Running time: %s Seconds'%(end-start))
       
    # remove the dud before saving
    dud_AAC = []
    for pdbid in contact:
        if contact[pdbid] == []:
            dud_AAC.append(pdbid)
    dud_AAC
    for dud in dud_AAC:
        del contact[dud]
        del combined_ids[dud]
        del matched_ids[dud]
           

        
    return sequence, contact, matched_ids, combined_ids, training_ids, testing_ids
 
    
   
######################################################################
if __name__ == '__main__':
    # wd is the working directory, sd is the saving directory
    wd = "/home/leo/Documents/Database/Pipeline_New/All with peptide 5+ resolution 4A/structure"
    sd = '/home/leo/Documents/Database/Pipeline_New/All with peptide 5+ resolution 4A'
    
    sequence, contact, matched_ids, combined_ids, training_ids, testing_ids =  main(wd, 0.1)
    
    # Separete the combined_ids and the matched_ids into the training and testing subsets
    training_combined_ids ={}
    training_matched_ids ={}
    testing_combined_ids = {}
    testing_matched_ids = {}
    
    for key, value in combined_ids.items():
        if key in training_ids:
            training_combined_ids[key] = value
        elif key in testing_ids:
            testing_combined_ids[key] = value
            
    for key, value in matched_ids.items():
        if key in training_ids:
            training_matched_ids[key] = value
        elif key in testing_ids:
            testing_matched_ids[key] = value
    
    # Separate the extracted sequences into training and testing subsets
    training_sequence = {}
    testing_sequence = {}
    for key, value in sequence.items():
        if key in training_ids:
            training_sequence[key] = value
        elif key in testing_ids:
            testing_sequence[key] = value
            
    # Separate the contact into training and testing subsets
    training_contact = {}
    testing_contact = {}
    for key, value in contact.items():
        if key in training_ids:
            training_contact[key] = value
        elif key in testing_ids:
            testing_contact[key] = value
            
    # Save the above results    
    os.chdir(sd)
    with open('sequence', 'w') as f:
        json.dump(sequence, f)
    with open('contact', 'w') as f:
        json.dump(contact, f)
    with open('matched_ids', 'w') as f:
        json.dump(matched_ids, f)
    with open('combined_ids', 'w') as f:
        json.dump(combined_ids, f)
    with open('training_ids', 'w') as f:
        json.dump(training_ids, f)
    with open('testing_ids', 'w') as f:
        json.dump(testing_ids, f)    
        
    with open('training_combined_ids', 'w') as f:
        json.dump(training_combined_ids, f)
    with open('testing_combined_ids', 'w') as f:
        json.dump(testing_combined_ids, f)
    with open('training_matched_ids', 'w') as f:
        json.dump(training_matched_ids, f)
    with open('testing_matched_ids', 'w') as f:
        json.dump(testing_matched_ids, f)
    with open('training_sequence', 'w') as f:
        json.dump(training_sequence, f)
    with open('testing_sequence', 'w') as f:
        json.dump(testing_sequence, f)
    with open('training_contact', 'w') as f:
        json.dump(training_contact, f)
    with open('testing_contact', 'w') as f:
        json.dump(testing_contact, f)

#len(sequence)       
#len(matched_ids)    
#len(combined_ids)
#len(training_ids)
#len(testing_ids)
#len(training_combined_ids)
#len(testing_combined_ids)
#len(training_sequence)
#len(testing_sequence)
#len(training_contact)
#len(testing_contact)
#len(training_matched_ids)
#len(testing_matched_ids)
#
#training_contact.keys()
#training_contact['1fsk']
