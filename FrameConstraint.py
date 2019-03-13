import copy
import os
import json
#########################################################
#os.chdir('/home/leo/Documents/Database/Pipeline_New/All with peptide 5+ resolution 4A')
#with open('training_ac_contact', 'r') as f:
#    ac_contact = json.load(f)
#with open('training_sequence', 'r') as f:
#    sequence = json.load(f)
#with open('contact', 'r') as f:
#    contact = json.load(f)


##################################################################################
'''
Add:
    a function to extend from the number header so a sequence, so that the 
    length of this sequence is len(interval)+1, and the difference between to 
    consecutive numbers equal the numbers given in the interval
     
    for example, Add([1,1,2], 2), the returned list should be
    [2, 3, 4, 6], where 3-2, 4-3, 6-4 are exactly the numbers in [1, 1, 2], and 
    the starting number is 2.
INput:
    interval:
        a list of positive integers
    hearder:
        an integer
Output:
     res:
         a list, as explained in the above.
'''
def Add(interval, header):
    res = [header]
    m = 0
    while len(res) <= len(interval):
        res.append(res[-1]+interval[m])
        m += 1
    return res
#############################################################################
'''
Sub_seq:
    a function to tell whether a sub_sequence is a subsequence of a given sequence
Input:
    sequence:
        a list of integers
    sub_sequence:
        a list of integers
Output:
    boolean:
        boolean, if the sub_sequence is a sub sequence of sequence, then it takes 
        the value of True. Otherwise, it takes the value of False.
'''

def Sub_seq(sequence, sub_sequence):
    boolean = True
    for i in sub_sequence:
        if i not in sequence:
            boolean = False
            break
    return boolean
#Sub_seq([1, 2, 3], [1, 4])
#####################################################################################
'''
Get_consecutive:
     a function to find all the consecutive numbers with given length and free_type
Input:
    sequence:
        a list of numbers 
    length:
        an integer, gives the length of the amino acids
    free_type:
        takes either 0 or 1, which means 0 free or 1 free.
Output:
    sub_sequence:
        a list, of which each element is a list of consecutive numbers with the given 
        length and free_type.
'''
def Get_consecutive(sequence, length, free_type):
    sequence.sort()
    sub_sequence = []
    if length == 1:
        for i in sequence:
            sub_sequence.append([i])
    else:
        # General the distance between two numbers, such that
        # the distance is no larger than the free_type
        interval = []
        for i in range(free_type+1):
            interval.append([i+1])
        while len(interval[0]) + 1 <length:
            longer_interval = []
            for inter in interval:            
                for j in range(free_type+1):
                    longer_inter = copy.deepcopy(inter)
                    longer_inter.append(j+1)
                    longer_interval.append(longer_inter)
            interval = longer_interval
        # Generate the subsequence with the given length and free_type 
        # We don't want to change the input sequence, thus we deep copy it.
        copy_sequence = copy.deepcopy(sequence)
        while len(copy_sequence) >= length:
            temp_sub_sequence = []
            for i in interval:
                temp_sub_sequence = Add(i, copy_sequence[0])
                if Sub_seq(copy_sequence, temp_sub_sequence):
                    if temp_sub_sequence not in sub_sequence:
                        sub_sequence.append(temp_sub_sequence)
            copy_sequence.remove(copy_sequence[0])        
    return sub_sequence

####################################################################################

'''
****** Explanations of the FC_parameter**************

    FC_parameter={}
    FC_parameter['ac_contact'] = ac_contact
    FC_parameter['contact'] = contact
    FC_parameter['sequence'] = sequence
    FC_parameter['CN_gate'] = CN_gate 
    FC_parameter['Ag_length'] = Ag_length
    FC_parameter['Ab_length'] = Ab_length
    FC_parameter['Ag_free_type'] = Ag_free_type
    FC_parameter['Ab_free_type'] = Ab_free_type
    FC_parameter['core_number_per_chain'] = core_number_per_chain
    FC_parameter['core_separation'] = core_separation
    FC_parameter['l_range'] = l_range  
    FC_parameter['h_range'] = h_range 

ac_contact:
    this is the output of the module AlignmentConstraint
sequence:
    this is the output of the module AAC_2
CN_gate:
    an integer, if the contact number of two amino acids is smaller than the CN_gate, they
    are viewed as 'no contact'
Ag_length:
    an integer, gives the length of the antigen amino acids
Ab_length:
    an integer, gives the length of the CDR amino acids length
Ag_free_type:
    an integer takes the value of either 0 or 1, which means either 0 free or 1 free.
    0_free, means 1,2 are next to each other but 1, 3 are not. 1 free means two integers
    are next to each other if the difference is no larger than 1.
core_number_per_chain:
    an integer, gives the number of cores per chain
core_separation:
    an integer, gives the number of amino acids two cores must be separated in order to be
    considered as two cores
l_range: 
    gives the locations of the CDRLs, refer to AAC_2
h_range:
    gives the locations of the CDRHs, refer to AAC_2
'''
###################################################################################

'''
CN_gated: 
    To get rid of the amino acid contact with contact number smaller than CN_gate
    
Output:
    parameter['cn_gated'], a dictionary in the same form as the ac_contact, the only difference
      is that the amino acids with contact number smaller than CN_gate is deleted
'''             
        
def CN_gated(FC_parameter):
    # Take out the value
    ac_contact = FC_parameter['ac_contact']
    CN_gate = FC_parameter['CN_gate']
    
    cn_gated = {}
    for Ab_chain in ac_contact:
        cn_gated[Ab_chain] = []
        for fcdn in ac_contact[Ab_chain]:
            if fcdn[3] >= CN_gate:
                cn_gated[Ab_chain].append(fcdn )
    FC_parameter['cn_gated'] = cn_gated  
#############################################################
'''
Ab_Ag_contact:
    a function to separate the cn_gated into groups according to different chain combination
'''
def Ab_Ag_Contact(FC_parameter):
    # Take out the values
    cn_gated = FC_parameter['cn_gated']
    
    Ab_Ag_contact = {}
    for key, value in cn_gated.items():
        for fcdn in value:
            if key+fcdn[0][3] not in Ab_Ag_contact:
                Ab_Ag_contact[key+fcdn[0][3]] = [fcdn]
            else:
                Ab_Ag_contact[key+fcdn[0][3]].append(fcdn)
    # load this to the FC_parameter
    FC_parameter['Ab_Ag_contact'] = Ab_Ag_contact
    
##############################################################
'''
Match_up_indices:
    a function to match up the consecutives amino acid positions of Ab with the
    consecutive amino acids positions of Ag, while all the consetive positions can
    meet the requirement
Output:
    FC_parameter['matched_up_indices'], a dictionary with the same keys as 
    FC_parameter['ac_contact']. the values are in the form of 
    [[[1, 2], [4, 6],'l2LA'], [[1,3], [6, 7], h3HA],...]
    where each element is in the form of [[1, 2], [4, 6]], which means the amino acids of CDR
    in position 1, 2 contact with amino acids of Ag in position 4, 6. And the positions meet the 
    requirements of length and free_type, 'l2LA' is the tag
'''
def Match_up_indices(FC_parameter):
    # Take out the values
    Ab_Ag_contact= FC_parameter['Ab_Ag_contact']
    Ab_length = FC_parameter['Ab_length']
    Ag_length = FC_parameter['Ag_length']
    Ab_free_type = FC_parameter['Ab_free_type']
    Ag_free_type = FC_parameter['Ag_free_type']
      
    # Create an empty dictionary to contain the results
    matched_up_indices = {}
    for key, value in Ab_Ag_contact.items():
        matched_up_indices[key] = []
        Ab_indx_sequence = []
        for fcdn in value:
            Ab_indx_sequence.append(fcdn[1])
            
        Ab_indx_sequence.sort()
        Ab_indx_sub_sequence = Get_consecutive(Ab_indx_sequence, Ab_length, Ab_free_type)
        
        # Match up the consecutive amino acids in for each consecutive amino acids in Ab
        if Ab_indx_sub_sequence != []:
            for Ab_sub in Ab_indx_sub_sequence:
                Ag_indx_sequence_for_Ab_sub = []
                for fcdn in value:
                    if fcdn[1] in Ab_sub:
                        Ag_indx_sequence_for_Ab_sub.append(fcdn[2])

                Ag_indx_sub_sequence_for_Ab_sub = Get_consecutive(\
                                                Ag_indx_sequence_for_Ab_sub, Ag_length, Ag_free_type)
            
                if Ag_indx_sub_sequence_for_Ab_sub != []:
                    for Ag_sub_for_Ab_sub in Ag_indx_sub_sequence_for_Ab_sub:
                        matched_up_indices[key].append([Ab_sub, Ag_sub_for_Ab_sub])
                
    FC_parameter['matched_up_indices'] = matched_up_indices
    

################################################################################## 
'''
Adjust_the_order:
    A function to keep the original order or reverse the orginal order, according 
    to the number of contact at the two ends of the matched up paires
Output:
    The order of the matched_up_indices are adjusted
'''

def Adjust_the_order(FC_parameter):
    # take out the values
    Ab_Ag_contact = FC_parameter['Ab_Ag_contact']
    matched_up_indices = FC_parameter['matched_up_indices']
    Ab_length = FC_parameter['Ab_length']
    Ag_length = FC_parameter['Ag_length']
    
    if Ab_length >= 2 and Ag_length >= 2:
    
        for key, value in matched_up_indices.items():
            for match in value:
                normal = 0
                reverse = 0
                for fcdn in Ab_Ag_contact[key]:
                    if fcdn[1] == match[0][0] and fcdn[2] == match[1][0] :
                        normal += fcdn[3]
                    elif fcdn[1] == match[0][-1] and fcdn[2] == match[1][-1]:
                        normal += fcdn[3]
                    elif fcdn[1] == match[0][0] and fcdn[2] == match[1][-1]:
                        reverse += fcdn[3]
                    elif fcdn[1] == match[0][-1] and fcdn[2] == match[1][0]:
                        reverse += fcdn[3]

                if reverse > normal:
                    match[0].sort(reverse=True)

    
################################################################################
'''
Core_selection:
    this function is to select the cores
Input:
    matched_up_indices_list:
        a list in the form of 
        [[[1, 2], [4, 5], 14, h3HA], ...]
    core_number_per_chain:
        an integer, give the maximum number of cores to be selected from the each Ab chain
    core_separation:
        the minimum distance to separate two cores
Output:
    selected_cores:
        a list in the same form as matched_up_indices_list
'''
def Core_selection(matched_up_indices_list, core_number_per_chain, core_separation):
    # Creat an empty container
    selected_cores = []
    # Sort the list according to the contact numbers
    matched_up_indices_list.sort(key = lambda x:x[2], reverse = True)
    # Make a deep copy
    value = copy.deepcopy(matched_up_indices_list)
    # Select
    while value!= [] and len(selected_cores) < core_number_per_chain:
        selected_cores.append(value[0])
        to_be_removed =[]
        for match in value:
            Ab_pos_start = min(match[0])
            Ab_pos_end = max(match[0])
            selected_Ab_pos_start = min(selected_cores[-1][0])
            selected_Ab_pos_end = max(selected_cores[-1][0])
#            selected_Ab_pos.sort()
            # make sure the distance between the selected Ab pos and the Ab pos in the value
            # are separated by at least core_separation number of amino acids
            if abs(Ab_pos_start-selected_Ab_pos_end)<= core_separation or \
            abs(Ab_pos_end - selected_Ab_pos_start) <= core_separation or \
            abs(Ab_pos_start - selected_Ab_pos_start) <= core_separation or\
            abs(Ab_pos_end - selected_Ab_pos_end) <= core_separation:
                to_be_removed.append(match)
        for i in to_be_removed:
            value.remove(i)
    return selected_cores
##################################################################################  

'''
Select_match_up_indices:
    a function to select the match_up with the largest contact and the order of the
    sequence adjusted.
Output:
    FC_parameter['core']:
        a list, with each element in the form of 
        [['ALA', 'THR'], ['PHE', 'GLU'], contact_number, l3LA, pdbid]
    FC_parameter['core_indx']:
        a list, with each element in the form of 
        [[15, 16], [37, 38], contact_number, l3, pdbid]
        
    ###The only difference between core and core_indx is that the triple letter
        is replaced by the position indices, for those position indices may be used
        for further analysis.###
'''
def Select_match_up_indices(FC_parameter):
    # Take out the values
#    cn_gated = FC_parameter['cn_gated']
    matched_up_indices = FC_parameter['matched_up_indices']
    Ab_Ag_contact = FC_parameter['Ab_Ag_contact']
    core_number_per_chain = FC_parameter['core_number_per_chain']
    core_separation = FC_parameter['core_separation']
    # Calculate the contact number for each combination, and tag along  with four letters 'l3LA'  
    # and slecte the core_number_per_chain of cores from each chain
    selected_match_up_indices = {}
    for key, value in matched_up_indices.items():
        selected_match_up_indices[key] = []
        if value!= []:
            for match in value:
                n_contact = 0
                cdr = ''
                for fcdn in Ab_Ag_contact[key]:
                    if fcdn[1] in match[0] and fcdn[2] in match[1]:
                        n_contact += fcdn[3]
                        cdr = fcdn[0]
                match.extend([n_contact, cdr, key[:4]])
                
            # Select the cores
            cores = Core_selection(value, core_number_per_chain, core_separation)
            selected_match_up_indices[key]=cores
        
    FC_parameter['selected_match_up_indices'] = selected_match_up_indices
########################################################################
'''
Choose_the_best:
    a small function to be used in the Each_CDR_get_one to choose the 
    match up with the largest contact number
'''
def Choose_the_best(G):
    max_contact = 0
    best = []
    for match in G:
        if len(match) < 3:
            print(match)
            break
        if match[2] > max_contact:
            max_contact = match[2]
            best = match
    return best                    
######################################################################
'''
Each_CDR_get_one:
    a function to try to select a core on each CDR
    
    ###This function can only be run after Select_match_up_indices, because
        it needs the 'matched_up_indices' be modified ###
Input: FC_parameter
Output:
    FC_parameter['each_CDR_get_one']:
        in the same form as the FC_parameter['selected_match_up_indices']
        
'''
def Each_CDR_get_one(FC_parameter):
    matched_up_indices = FC_parameter['matched_up_indices']
    l_range = FC_parameter['l_range']
    h_range = FC_parameter['h_range']
    
    each_CDR_get_one = {}
    for key, value in matched_up_indices.items():
        each_CDR_get_one[key] = []
        if key[5] == 'l':
            Rs = l_range
        elif key[5] == 'h':
            Rs = h_range
        R1=Rs[0];R2=Rs[1];R3=Rs[2]  
        # Group the value into three groups according to the Rs
        G1=[]; G2=[]; G3=[]
        for match in value:
            if match[0][0]>= R1[0] and match[0][0] <= R1[-1]:
                G1.append(match)
            elif match[0][0]>= R2[0] and match[0][0] <= R2[-1]:
                G2.append(match)
            elif match[0][0]>= R3[0] and match[0][0] <= R3[-1]:
                G3.append(match)
        # Select the one with the largest contact number from each group
        each_CDR_get_one[key] = [Choose_the_best(G1), Choose_the_best(G2), Choose_the_best(G3)]
        
    FC_parameter['each_CDR_get_one'] = each_CDR_get_one

########################################################################
'''
Change_to_aa:
    This function is to find the coresponding amino acids according to the information
    of the indices given in the FC_parameter['selected_match_up_indices'].
'''

def Change_to_aa(FC_parameter):
    # take out the values
    selected_match_up_indices = FC_parameter['selected_match_up_indices']
    sequence = FC_parameter['sequence']
    
    # Change to aa
    core_aa = []
    for key, value in selected_match_up_indices.items():
        chain_sequence = sequence[key[:4]]
        for match in value:
            if match != []:
                Ab_pos = match[0]
                Ag_pos = match[1]
                n_contact = match[2]
                cdl = match[3]
                
                Ab_aa = []            
                for i in Ab_pos:
                    Ab_aa.append(chain_sequence[cdl[2]][i])
                    
                Ag_aa = []
                for j in Ag_pos:
                    Ag_aa.append(chain_sequence[cdl[3]][j])
                
                core_aa.append([Ab_aa, Ag_aa, n_contact, cdl, match[4]])
            
    FC_parameter['core_aa'] = core_aa
        
  
######################################################
'''
The_middle_aa:
    This function is to find the middle aa if the selected_match_up_indices are not 
    next to each other. This function is trivial when the free_type is 0.
'''
def The_middle_aa(FC_parameter):
    # Take out the values
    Ab_Ag_contact = FC_parameter['Ab_Ag_contact']
    sequence = FC_parameter['sequence']
    
    # Change to aa
    middle_Ab_aa = []
    middle_Ag_aa = []
    for key, value in Ab_Ag_contact.items():
        chain_sequence = sequence[key[:4]]
        Ab_pos = []
        Ag_pos = []
        for fcdn in value:
            Ab_pos.append(fcdn[1])
            Ag_pos.append(fcdn[2])
            
        Ab_pos.sort()
        Ag_pos.sort()
                    
        for i in range(len(Ab_pos)-1):
            if Ab_pos[i+1]-Ab_pos[i] == 2:
                flank_middle = chain_sequence[key[4]][Ab_pos[i]:Ab_pos[i]+3]
                flank_middle.append(key)
                middle_Ab_aa.append(flank_middle)
        
        for i in range(len(Ag_pos)-1):
            if Ag_pos[i+1]-Ag_pos[i] == 2:
                flank_middle = chain_sequence[key[6]][Ag_pos[i]:Ag_pos[i]+3]
                flank_middle.append(key)
                middle_Ag_aa.append(flank_middle)

           
    FC_parameter['middle_Ab_aa'] = middle_Ab_aa
    FC_parameter['middle_Ag_aa'] = middle_Ag_aa
         
################################################################################            
           
def Main():
    
    # Set the working directory and the saving directory
    wd = '/home/leo/Documents/Database/Pipeline_New/All with peptide 5+ resolution 4A'
    sd = '/home/leo/Documents/Database/Pipeline_New/Cores'

    
    training_testing = ['training', 'testing']
    for train_test in training_testing:
        # Go to the working directory
        os.chdir(wd)
        with open(train_test+'_ac_contact', 'r') as f:
            ac_contact = json.load(f)
        with open(train_test+'_sequence', 'r') as f:
            sequence = json.load(f)
        with open('contact', 'r') as f:
            contact = json.load(f)
            
        FC_parameter={}
        FC_parameter['ac_contact'] = ac_contact
        FC_parameter['contact'] = contact
        FC_parameter['sequence'] = sequence   
         
        CN_gate = 1
        FC_parameter['CN_gate'] = CN_gate 
        
        CN_gated(FC_parameter) 
        Ab_Ag_Contact(FC_parameter)
        # Calculate the middle aa
        The_middle_aa(FC_parameter)
        # Save the results
        with open(train_test+'_middle_Ab_aa', 'w') as f:
            json.dump(FC_parameter['middle_Ab_aa'], f)
        with open(train_test+'_middle_Ag_aa', 'w') as f:
            json.dump(FC_parameter['middle_Ag_aa'], f)
        # Select the cores by assigning different number of cores per chain
        for i in range(1, 5):
            for j in range(1, 5):                
                Ag_length =i
                Ab_length = j
                FC_parameter['Ag_length'] = Ag_length
                FC_parameter['Ab_length'] = Ab_length
                
                for k in [1, 0]:
                    Ag_free_type = k
                    Ab_free_type = k
                    FC_parameter['Ag_free_type'] = Ag_free_type
                    FC_parameter['Ab_free_type'] = Ab_free_type
                    
                    l_range = [[23, 40], [49, 63], [89, 110]]
                    h_range = [[25, 37], [50, 71], [99, 129]]
                    FC_parameter['l_range'] = l_range  
                    FC_parameter['h_range'] = h_range 
                    
                    Match_up_indices(FC_parameter)
                    Adjust_the_order(FC_parameter)
                        
                    for h in [1, 2, 3]:
                        core_number_per_chain = h
                        core_separation = 2                       
                        FC_parameter['core_number_per_chain'] = core_number_per_chain
                        FC_parameter['core_separation'] = core_separation
                        
                        Select_match_up_indices(FC_parameter)
                        Change_to_aa(FC_parameter)
                                                                   
                        # Save the results
#                        os.chdir(wd)
                        name = str(i)+'_'+str(j)+'_'+str(k)+'_'+str(k)+'_1_2_'+str(h)+'perchain'
                        print('working on '+ name)

                        # Save the core aa
                        os.chdir(sd)
                        with open(train_test+'_'+name, 'w') as f:
                            json.dump(FC_parameter['core_aa'], f)
                            
                    Each_CDR_get_one(FC_parameter)
                    FC_parameter['selected_match_up_indices'] = FC_parameter['each_CDR_get_one']
                    Change_to_aa(FC_parameter)                    
                    # save the results
                    os.chdir(sd)
                    name = str(i)+'_'+str(j)+'_'+str(k)+'_'+str(k)+'_1_'+'1perCDR'
                    print('working on '+ name)
                    with open(train_test+'_'+name, 'w') as f:
                        json.dump(FC_parameter['core_aa'], f)

    
##########################################################################
#if __name__ == '__main__':
#    Main()

                            
#############################################################################
'''
# Set the working directory and the saving directory
wd = '/home/leo/Documents/Database/Pipeline_New/All with peptide 5+ resolution 4A'
sd = '/home/leo/Documents/Database/Pipeline_New/Cores'


training_testing = ['training', 'testing']
for train_test in training_testing:
    # Go to the working directory
    os.chdir(wd)
    with open(train_test+'_ac_contact', 'r') as f:
        ac_contact = json.load(f)
    with open(train_test+'_sequence', 'r') as f:
        sequence = json.load(f)
    with open('contact', 'r') as f:
        contact = json.load(f)
        
    FC_parameter={}
    FC_parameter['ac_contact'] = ac_contact
    FC_parameter['contact'] = contact
    FC_parameter['sequence'] = sequence   
     
    CN_gate = 1
    FC_parameter['CN_gate'] = CN_gate 
    
    CN_gated(FC_parameter) 
    Ab_Ag_Contact(FC_parameter)
    # Calculate the middle aa
#    The_middle_aa(FC_parameter)
#    # Save the results
#    with open(train_test+'_middle_Ab_aa', 'w') as f:
#        json.dump(FC_parameter['middle_Ab_aa'], f)
#    with open(train_test+'_middle_Ag_aa', 'w') as f:
#        json.dump(FC_parameter['middle_Ag_aa'], f)
    # Select the cores by assigning different number of cores per chain
    for i in range(1, 5):
        for j in range(1, 5):                
            Ag_length =i
            Ab_length = j
            FC_parameter['Ag_length'] = Ag_length
            FC_parameter['Ab_length'] = Ab_length
            
            for k in [1, 0]:
                Ag_free_type = k
                Ab_free_type = k
                FC_parameter['Ag_free_type'] = Ag_free_type
                FC_parameter['Ab_free_type'] = Ab_free_type
                
                l_range = [[23, 40], [49, 63], [89, 110]]
                h_range = [[25, 37], [50, 71], [99, 129]]
                FC_parameter['l_range'] = l_range  
                FC_parameter['h_range'] = h_range 
                
                Match_up_indices(FC_parameter)
                Adjust_the_order(FC_parameter)
                    
                for h in [1, 2, 3]:
                    core_number_per_chain = h
                    core_separation = 2                       
                    FC_parameter['core_number_per_chain'] = core_number_per_chain
                    FC_parameter['core_separation'] = core_separation
                    
                    Select_match_up_indices(FC_parameter)
                    Change_to_aa(FC_parameter)
                                                               
                    # Save the results
#                        os.chdir(wd)
                    name = str(i)+'_'+str(j)+'_'+str(k)+'_'+str(k)+'_1_2_'+str(h)+'perchain'
#                    print('working on '+ name)
                    

                    print (train_test+'_'+name, len(FC_parameter['core_aa']))

#                    # Save the core aa
#                    os.chdir(sd)
#                    with open(train_test+'_'+name, 'w') as f:
#                        json.dump(FC_parameter['core_aa'], f)
                        
                Each_CDR_get_one(FC_parameter)
                FC_parameter['selected_match_up_indices'] = FC_parameter['each_CDR_get_one']
                Change_to_aa(FC_parameter)                    
                # save the results
                os.chdir(sd)
                name = str(i)+'_'+str(j)+'_'+str(k)+'_'+str(k)+'_1_'+'1perCDR'
#                print('working on '+ name)
                print(train_test+'_'+name, len(FC_parameter['core_aa'])) 

#                with open(train_test+'_'+name, 'w') as f:
#                    json.dump(FC_parameter['core_aa'], f)
############################################################################
#'''
#The following comments is for testing the above function with a small artifical data set.
#'''
#
#
## Create a small testing data
#test_sequence = {}
#test_sequence['test'] = {}
#test_sequence['test']['H'] = ['SER','ASP','SER','LYS','ILE','LEU','ALA','HIS','LEU','PHE', 'ILE', 'SER', 'LYS']
#test_sequence['test']['L'] = ['GLY','TYR','ASP','PHE','ARG','VAL','ARG','PRO','PRO','THR', 'ILE','LEU','ALA']
#test_sequence['test']['A'] = ['ASP','ASN','GLY','GLY','PRO','ASN','MET','LEU','LEU','ARG','PHE','ARG','VAL']
#test_sequence['test']['B'] = ['ASP','ASN','GLY','GLY','PRO','ASN','MET','LEU','LEU','ARG','PHE','ARG','VAL']
#
#test_contact = {}
#test_contact['test'] = [['l1LA', 0, 2, 5],
#                            ['l1LA', 2, 1, 8],
#                            ['l2LA', 5, 4, 3],
#                            ['l2LA', 6, 6, 2],
#                            ['l3LA', 10, 7,8],
#                            ['l3LA', 11, 8, 8],
#                            ['h1HA', 0, 0, 1], 
#                            ['h1HA', 1, 0, 3],
#                            ['h2HA', 5, 8, 4],
#                            ['h2HA', 6, 7, 5],
#                            ['h3HA', 10, 5, 10],
#                            ['h3HA', 12, 7, 9],
#                            ['h1HB', 0, 0, 1], 
#                            ['h1HB', 1, 0, 3],
#                            ['h2HB', 5, 8, 4],
#                            ['h2HB', 6, 7, 5],
#                            ['h3HB', 10, 5, 10],
#                            ['h3HB', 12, 7, 9]]
#test_ac_contact={}
#test_ac_contact['testHh']=  [['h1HA', 0, 0, 1], 
#                            ['h1HA', 1, 0, 3],
#                            ['h2HA', 5, 8, 4],
#                            ['h2HA', 6, 7, 5],
#                            ['h3HA', 10, 5, 10],
#                            ['h3HA', 12, 7, 9],
#                            ['h1HB', 0, 0, 1], 
#                            ['h1HB', 1, 0, 3],
#                            ['h2HB', 5, 8, 4],
#                            ['h2HB', 6, 7, 5],
#                            ['h3HB', 10, 5, 10],
#                            ['h3HB', 12, 7, 9]]
#
#test_ac_contact['testLl'] =  [['l1LA', 0, 2, 5],
#                            ['l1LA', 2, 1, 8],
#                            ['l2LA', 5, 4, 3],
#                            ['l2LA', 6, 6, 2],
#                            ['l3LA', 10, 7,8],
#                            ['l3LA', 11, 8, 8],
#                            ['l3LA', 12, 6, 9]]
#
#             
#
#CN_gate = 1
#Ag_length =2
#Ab_length = 2
#Ag_free_type = 0
#Ab_free_type = 0
#core_number_per_chain = 1
#core_separation = 2
#l_range = [[0, 2], [5,7], [10, 12]]
#h_range = [[0, 2], [5,7], [10, 12]]
#FC_parameter={}
#FC_parameter['ac_contact'] = test_ac_contact
#FC_parameter['contact'] = test_contact
#FC_parameter['sequence'] = test_sequence
#FC_parameter['CN_gate'] = CN_gate 
#FC_parameter['Ag_length'] = Ag_length
#FC_parameter['Ab_length'] = Ab_length
#FC_parameter['Ag_free_type'] = Ag_free_type
#FC_parameter['Ab_free_type'] = Ab_free_type
#FC_parameter['core_number_per_chain'] = core_number_per_chain
#FC_parameter['core_separation'] = core_separation
#FC_parameter['l_range'] = l_range  
#FC_parameter['h_range'] = h_range  
#
## Test the CN-gated
#FC_parameter['CN_gate'] = 1
#CN_gated(FC_parameter) 
#FC_parameter['cn_gated']
##Test the Match_up_indices
#FC_parameter['Ag_length'] = Ag_length
#FC_parameter['Ab_length'] = Ab_length
#FC_parameter['Ag_free_type'] = 1
#FC_parameter['Ab_free_type'] = 1
#Match_up_indices(FC_parameter)
#FC_parameter['matched_up_indices']
## Test Adjust_the_order
#Adjust_the_order(FC_parameter) 
#FC_parameter['matched_up_indices']
## Test Select_match_up_indices
#FC_parameter['core_number_per_chain'] = 2
#FC_parameter['core_separation'] = 2
#Select_match_up_indices(FC_parameter)
##FC_parameter.keys()
#FC_parameter['selected_match_up_indices']
#FC_parameter['matched_up_indices']
## Test the Change_to_aa
#Change_to_aa(FC_parameter)
#FC_parameter['core_aa']
#The_middle_aa(FC_parameter)
#FC_parameter['middle_Ab_aa']
#FC_parameter['middle_Ag_aa']
#Each_CDR_get_one(FC_parameter)
#FC_parameter['each_CDR_get_one']
#FC_parameter['selected_match_up_indices'] = FC_parameter['each_CDR_get_one']
#Change_to_aa(FC_parameter)
#FC_parameter['core_aa']
  
