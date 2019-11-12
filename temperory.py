import json
import os
import math
from scipy.stats import norm



####################################################################################3
# Take a look at the distribution of cores over different CDRs
def Count_cores_over_CDRs():
    l_range = [[23, 40], [49, 63], [89, 110]]
    h_range = [[25, 37], [50, 71], [99, 129]]

    l_total = 0; h_total = 0
    for l, h in zip(l_range, h_range):
        l_total += l[1] - l[0] +1
        h_total += h[1] - h[0] +1
        
    l_prop = []; h_prop = []    
    for l, h in zip(l_range, h_range):
        l_prop.append((l[1]-l[0]+1)/l_total)
        h_prop.append((h[1]-h[0]+1)/h_total)
        
    
    l_exp_ob = {}; h_exp_ob = {}
    os.chdir('/home/leo/Documents/Database/Data_Code_Publish/Cores/Positive_cores')
    core_count_CDRs = {}
    for i in range(1,4):
        for j in range(1, 4):
            match_type = str(i)+'_'+str(j)
    
            test_name = 'testing_'+match_type+'_0_0_1_2_1perchain'
            train_name = 'training_'+match_type+'_0_0_1_2_1perchain'
            
            with open(test_name,'r') as f:
                testing = json.load(f)
            with open(train_name,'r') as f:
                training= json.load(f)
        
            cores_match_type = testing[:]
            cores_match_type.extend(training)
            
            h1 = 0; h2 = 0; h3 = 0; l1 = 0; l2 = 0; l3 = 0
            for fcdn in cores_match_type:
                if fcdn[3][:2] == 'h1':
                    h1 += 1
                elif fcdn[3][:2] == 'h2':
                    h2 += 1 
                elif fcdn[3][:2] == 'h3':
                    h3 += 1 
                elif fcdn[3][:2] == 'l1':
                    l1 += 1 
                elif fcdn[3][:2] == 'l2':
                    l2 += 1 
                elif fcdn[3][:2] == 'l3':
                    l3 += 1 
                    
            core_count_CDRs[match_type] = [(h1, 'h1'),(h2, 'h2'),(h3, 'h3'),(l1, 'l1'),\
                           (l2, 'l2'),(l3, 'l3')]
            
            l_exp_ob[match_type] = []
            h_exp_ob[match_type] = []
            
            l_total = l1 + l2 + l3
            l_exp_ob[match_type].append([l_total*l_prop[0], l_total*l_prop[1], l_total*l_prop[2]])
            l_exp_ob[match_type].append([l1, l2, l3])
            
            h_total = h1 + h2 + h3
            h_exp_ob[match_type].append([h_total*h_prop[0], h_total*h_prop[1], h_total*h_prop[2]])
            h_exp_ob[match_type].append([h1, h2, h3])

            
    return core_count_CDRs, l_exp_ob, h_exp_ob, l_prop, h_prop

#core_count_CDRs, l_exp_ob, h_exp_ob, l_prop, h_prop = Count_cores_over_CDRs()
#l_exp_ob.keys()
#l_prop
#h_prop
# Calculate the p values to tell whether proportion of cores in CDR3 are larger
# after corrected by the length
def One_p_z_test(l_exp_ob, l_prop):
    exp_prop = l_prop[2]
    pvs = {}
    for i in range(1,4):
        for j in range(1, 4):
            match_type = str(i)+'_'+str(j)
            ob = l_exp_ob[match_type][1]
            cdr_total = ob[0]+ob[1]+ob[2]
            ob_prop = ob[2]/cdr_total
            z = ob_prop-exp_prop
            z /= (exp_prop*(1-exp_prop)/cdr_total)**0.5
            p_value = 1 - norm.cdf(z)
            
            pvs[match_type] = p_value
    return pvs


# Pack up the above functions
def Core_over_CDR_statistics():
    core_count_CDRs, l_exp_ob, h_exp_ob, l_prop, h_prop = Count_cores_over_CDRs()
    chains = ['heavy_chain', 'light_chain']
    p_values_CDR3_more_cores = {}
    for chain in chains:
        if chain == 'heavy_chain':
            p_values_CDR3_more_cores[chain]=One_p_z_test(h_exp_ob, h_prop)
        elif chain == 'light_chain':
            p_values_CDR3_more_cores[chain] = One_p_z_test(l_exp_ob, l_prop)
    # Load all the results into a dicitonary
    core_over_CDR_statistics = {}
    core_over_CDR_statistics['p_values_CDR3_more_cores'] = p_values_CDR3_more_cores
    core_over_CDR_statistics['core_count_CDRs'] = core_count_CDRs            
    
    return core_over_CDR_statistics

'''##########################################################################'''
# Calculate the statistics of the distribution of amino acids
os.chdir('/home/leo/Documents/Database/Data_Code_Publish/Structures')
with open('sequence', 'r') as f:
    sequences = json.load(f)
with open('combined_ids', 'r') as f:
    combined_ids = json.load(f)
with open('matched_ids', 'r') as f:
    matched_ids = json.load(f)

def Top_5_percent_cores():
    os.chdir('/home/leo/Documents/Database/Data_Code_Publish/Cores/Positive_cores')
    top_5_percent_cores = {}
    for i in range(1, 4):
        for j in range(1,4):
            match_type = str(i)+'_'+str(j)
            test_name = 'testing_'+match_type+'_0_0_1_2_1perchain'
            train_name = 'training_'+match_type+'_0_0_1_2_1perchain'
            
            
            with open(test_name, 'r') as f:
                test = json.load(f)
            with open(train_name, 'r') as f:
                train = json.load(f)
            
            # Extract the top 5 percent cores, frequency and relative frequency
            cores_list = test[:]
            cores_list.extend(train)
            
            cores = []
            for fcdn in cores_list:
                AA=''
                for Ab_aa in fcdn[0]:
                    AA += Ab_aa
                for Ag_aa in fcdn[1]:
                    AA += Ag_aa
                cores.append(AA)
                
            cores_set = set(cores)
            
            cores_count = []
            for c in cores_set:
                cores_count.append([c, cores.count(c)])
            cores_count.sort(key=lambda x:x[1], reverse=True)
            
            # Append the relative frequency
            total = len(cores_list)
            for cc in cores_count:
                cc.append(round(cc[1]/total, 5))
            # Load to the dictionary
            top_5_percent_cores[match_type]=cores_count[:math.ceil(0.05*len(cores_set))]
                
    return top_5_percent_cores

# 
top_5_percent_cores = Top_5_percent_cores()





































    
