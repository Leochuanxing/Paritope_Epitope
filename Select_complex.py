'''
This file is to select the complexes with given resolution
'''
import os
import json
import copy
import math
'''
Resolution_restricted: a function to select the complexes with a given resolution or better

INput: 
    summary: a summary file in the form of tsv
Output:
    complex_list: a list of complexes with the id pdbid and the resolution the given resolution 
    or better
'''
def Resolution_restricted(resolution = 2.5):  
    os.chdir('/home/leo/Documents/Database/Pipeline_New/All with peptide 5+ resolution 4A/structure')
    file = open('summary.tsv', 'r')
    summary = file.readlines()
    file.close()
    # Creat a container
    complex_list = []
    for i in range(1, len(summary)):
        line = summary[i].split('\t')
        if len(line) > 16 and float(line[16]) <= resolution:
            if line[0] not in complex_list:
                complex_list.append(line[0])
    
    return complex_list
#complex_list = Resolution_restricted(resolution=2.5)
#len(complex_list)

# Processing the latest data
def Resolution_restricted_latest(resolution = 2.5):
    os.chdir('/home/leo/Documents/Database/Pipeline_New/Latest')
    file_latest = open('summary.tsv', 'r')
    summary_latest = file_latest.readlines()
    file_latest.close()
    with open('combined_ids_latest', 'r') as f:
        good_combined_ids_latest = json.load(f)
    
    latest_ids = list(good_combined_ids_latest.keys())
    # Creat a container
    complex_list_latest = []
    for i in range(1, len(summary_latest)):
        line = summary_latest[i].split('\t')
        if len(line) > 16 and line[0] in latest_ids and float(line[16]) <= resolution:
            if line[0] not in complex_list_latest:
                complex_list_latest.append(line[0]) 
        
    return complex_list_latest

#complex_list_latest = Resolution_restricted_latest(resolution = 2.5)
#len(complex_list_latest)

def Number_date(date):
    date_number = ''
    date_list = date.split('/')
    if int(date_list[2]) < 20:
        date_number += '20'
    else: date_number += '19'
    date_number += date_list[2]
    date_number += date_list[0]
    date_number += date_list[1]
        
    date_number = int(date_number)
    
    return date_number



def Attach_date_to_complex(complex_list, complex_list_latest):
    os.chdir('/home/leo/Documents/Database/Pipeline_New/Latest')
    file_latest = open('summary.tsv', 'r')
    summary_latest = file_latest.readlines()
    file_latest.close()
    
    os.chdir('/home/leo/Documents/Database/Pipeline_New/All with peptide 5+ resolution 4A/structure')
    file = open('summary.tsv', 'r')
    summary = file.readlines()
    file.close()
    
    comp = copy.deepcopy(complex_list)
    comp.extend(complex_list_latest)
    
    comp_date = []
    for i in range(1, len(summary_latest)):
        line = summary_latest[i].split('\t')
        if line[0] in comp and [line[0], Number_date(line[9])] not in comp_date:
            comp_date.append([line[0], Number_date(line[9])])
            
    for j in range(1, len(summary)):
        line = summary[j].split('\t')
        if line[0] in comp and [line[0], Number_date(line[9])] not in comp_date:
            comp_date.append([line[0], Number_date(line[9])])
            
    return comp_date

#comp_date = Attach_date_to_complex(complex_list, complex_list_latest)


def Get_training_testing_complex(resolution=2.5, percentage=0.1):
    complex_list = Resolution_restricted(resolution)
    complex_list_latest = Resolution_restricted_latest(resolution)
    comp_date = Attach_date_to_complex(complex_list, complex_list_latest)
    
    # Put the latest percentage away as the testing set
    comp_date.sort(key = lambda x:x[1], reverse= True)
    n = math.floor(len(comp_date)*0.1)
    testing_set = comp_date[:n]
    training_set = comp_date[n:]
    
    return testing_set, training_set

   
def Get_combined_matched_ids(testing_set, training_set):
        
    os.chdir('/home/leo/Documents/Database/Pipeline_New/Latest')
    with open('combined_ids_latest', 'r') as f:
        combined_ids_latest = json.load(f)
    with open('matched_ids_latest', 'r') as f:
        matched_ids_latest = json.load(f)
        
    os.chdir('/home/leo/Documents/Database/Pipeline_New/All with peptide 5+ resolution 4A')
    with open('combined_ids', 'r') as f:
        combined_ids = json.load(f)
    with open('matched_ids', 'r') as f:
        matched_ids = json.load(f)
    
    for key, value in combined_ids_latest.items():
        combined_ids[key] =  value
    for k, v in matched_ids_latest.items():
        matched_ids[k] = v
    
    testing_combined_ids = {}; testing_matched_ids = {}
    training_combined_ids = {}; training_matched_ids = {}
    
    for com in testing_set:
        if com[0] in combined_ids:
            testing_combined_ids[com[0]] = combined_ids[com[0]]
        if com[0] in matched_ids:
            testing_matched_ids[com[0]] = matched_ids[com[0]]
    
    for comp in training_set:
        if comp[0] in combined_ids:
            training_combined_ids[comp[0]] = combined_ids[comp[0]]
        if comp[0] in matched_ids:
            training_matched_ids[comp[0]] = matched_ids[comp[0]]
            
    return testing_combined_ids, testing_matched_ids, training_combined_ids, training_matched_ids

def main(resolution, percentage):
    testing_set, training_set = Get_training_testing_complex(resolution, percentage) 
    testing_combined_ids, testing_matched_ids,\
    training_combined_ids, training_matched_ids = Get_combined_matched_ids(testing_set, training_set)
    
    return testing_set, training_set,\
             testing_combined_ids, testing_matched_ids,\
             training_combined_ids, training_matched_ids

if __name__ == '__main__':    
    testing_set, training_set,\
    testing_combined_ids, testing_matched_ids,\
    training_combined_ids, training_matched_ids = main(resolution=2.5, percentage=0.1)
    # Save the results
    os.chdir('/home/leo/Documents/Database/Pipeline_New/Complexes')
    with open('testing_set', 'w') as f:
        json.dump(testing_set, f)
    with open('training_set' , 'w') as f:
        json.dump(training_set, f)
    with open('testing_combined_ids', 'w') as f:
        json.dump(testing_combined_ids, f)
    with open('testing_matched_ids', 'w') as f:
        json.dump(testing_matched_ids, f)
    with open('training_combined_ids', 'w') as f:
        json.dump(training_combined_ids, f)
    with open('training_matched_ids', 'w') as f:
        json.dump(training_matched_ids, f)
        
    combined_ids = {};  matched_ids = {}
    for k, v in testing_combined_ids.items():
        combined_ids[k] = v
    for k, v in training_combined_ids.items():
        combined_ids[k] = v
    
    for k, v in testing_matched_ids.items():
        matched_ids[k] = v
    for k, v in training_matched_ids.items():
        matched_ids[k] = v
    with open('combined_ids', 'w') as f:
        json.dump(combined_ids, f)
    with open('matched_ids', 'w') as f:
        json.dump(matched_ids, f)
        
    os.chdir('/home/leo/Documents/Database/Pipeline_New/All with peptide 5+ resolution 4A')
    with open('sequence', 'r') as f:
        sequence = json.load(f)
    with open('contact', 'r') as f:
        contact = json.load(f)
    os.chdir('/home/leo/Documents/Database/Pipeline_New/Latest')
    with open('sequence', 'r') as f:
        sequence_latest = json.load(f)
    with open('contact_latest', 'r') as f:
        contact_latest = json.load(f)
        
    sequence.update(sequence_latest)

    contact.update(contact_latest)
    os.chdir('/home/leo/Documents/Database/Pipeline_New/Complexes')
    with open('sequence', 'w') as f:
        json.dump(sequence, f)
    with open('contact', 'w') as f:
        json.dump(contact, f)
        
    
#os.chdir('/home/leo/Documents/Database/Pipeline_New/Complexes')
#with open('testing_set', 'r') as f:
#    testing_set  = json.load(f)
#with open('training_set' , 'r') as f:
#    training_set = json.load(f)
#with open('testing_combined_ids', 'r') as f:
#    testing_combined_ids = json.load(f)
#with open('testing_matched_ids', 'r') as f:
#    testing_matched_ids  = json.load(f)
#with open('training_combined_ids', 'r') as f:
#    training_combined_ids = json.load(f)
#with open('training_matched_ids', 'r') as f:
#    training_matched_ids = json.load(f)
#with open('combined_ids', 'r') as f:
#    combined_ids = json.load(f)
#with open('matched_ids', 'r') as f:
#    matched_ids = json.load(f)
#with open('sequence', 'r') as f:
#    sequence = json.load(f)
#with open('contact', 'r') as f:
#    contact = json.load(f)
#
#    
#len(testing_set)
#len(testing_combined_ids)
#len(testing_matched_ids)
#len(training_set)
#len(training_combined_ids)
#len(training_matched_ids)
#len(combined_ids)
#len(matched_ids)
#len(contact)
#len(sequence)










