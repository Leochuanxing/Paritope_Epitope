import json
import os
from scipy.stats import norm


# Take a look at the dates
#os.chdir('/home/leo/Documents/Database/Data_Code_Publish/Structures')
#with open('training_combined_ids', 'r') as f:
#    training_combined_ids = json.load(f)
#with open('testing_combined_ids', 'r') as f:
#    testing_combined_ids = json.load(f)
#with open('train_test_ids_dates', 'r') as f:
#    train_test_ids_dates = json.load(f)
#    
#len(train_test_ids_dates)
#len(training_combined_ids)
#len(testing_combined_ids)
#train_test_ids_dates.keys()
#len(train_test_ids_dates['training_ids_dates'])
#len(train_test_ids_dates['testing_ids_dates'])
#train_test_ids_dates['training_ids_dates'][-5:]
#train_test_ids_dates['testing_ids_dates'][:5]
#train_test_ids_dates['testing_ids_dates'][-5:]
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












    
