import os
import json
import copy
##########################################################################
'''
To_seq:
    a functionto change amino acids sequences in triple letter notation into seq objects, which
    are used in alignment
Input:
    aa_sequence:
        a list of triple letter notaions such as ['TRP', 'SER'....]
Output:
    seq_obj:
        a seq object corresponding to the aa_sequence
'''

def To_seq(aa_sequence):
    from Bio.Seq import Seq
    from Bio.Alphabet import IUPAC
    # TripleSingle is a list of correspondence between the triple letter notation and the single letter notation
    TripleSingle =  [['TYR', 'Y'], ['LYS', 'K'],['ASP', 'D'], ['ASN', 'N'], ['TRP', 'W'], ['PHE', 'F'], ['GLN', 'Q'],
       ['GLU', 'E'], ['PRO', 'P'], ['GLY', 'G'], ['THR', 'T'],['SER', 'S'], ['ARG', 'R'], ['HIS', 'H'],
       ['LEU', 'L'], ['ILE', 'I'], ['CYS', 'C'], ['ALA', 'A'], ['MET', 'M'], ['VAL', 'V']]
    seq_obj = None

    seq_single_letter = ''
    for aa in aa_sequence:
        for TS in TripleSingle:
            if TS[0] == aa:
                seq_single_letter += TS[1]
    seq_obj = Seq(seq_single_letter, IUPAC.protein)
    
    return seq_obj
'''
SeqCDR:
    to convert the sequence into sequence object and be prepared for futher alignment
Input:
    sequence:
        a dictionary, which contains extracted sequences from the AAC_2. The sequence output of
        AAC_2 can be used here directly.
    matched_ids:
        in the same form with the output of AAC_2.
        
        ###The pdbids of the sequence and the matche_ids should be corresponding with each other###
    
Output:
    in the form of a list [['1abvhH', Seq_obj], ['1abvlL', Seq_obj]]
    ['1abvhH', Seq_obj] means the pdbid is 1abv heavy chain with chain name H, and 
    the Seq_obj is the sequence object of the concatenated CDR1, CDR2, CDR3 of chain H.
'''   
def SeqCDR(sequence, matched_ids):
    # Concatenate the sequences 
    seqCDRH = []
    seqCDRL = []
    for pdbid in matched_ids:
        for combination in matched_ids[pdbid]:
            if combination[2] != '':
                if combination[0] != '':
                    # Calculate the total length of the heavy chain
                    total_length_H = len(sequence[pdbid][combination[0]])
                    # Find the smaller value between the total_length and 129
                    end_CDRH3 = min(total_length_H, 130)
                    CDRH = [pdbid+'h'+combination[0]]
                    CDRH.append(To_seq(sequence[pdbid][combination[0]][25:38]))
                    CDRH.append(To_seq(sequence[pdbid][combination[0]][50:72]))
                    CDRH.append(To_seq(sequence[pdbid][combination[0]][99:end_CDRH3]))
                    seqCDRH.append(CDRH)
                if combination[1] != '':
                    # Calculate the total length of the light chain
                    total_length_L = len(sequence[pdbid][combination[1]])
                    # Find the smaller value between the total_length_L- and 107
                    end_CDRL3 = min(total_length_L, 108)
                    CDRL = [pdbid+'l'+combination[1]]
                    CDRL.append(To_seq(sequence[pdbid][combination[1]][23:40]))
                    CDRL.append(To_seq(sequence[pdbid][combination[1]][49:63]))
                    CDRL.append(To_seq(sequence[pdbid][combination[1]][89:end_CDRL3]))
                    seqCDRL.append(CDRL)                        

    return seqCDRH, seqCDRL

    
'''
Hcluster:
    A function to do hierarchical clustering of seqCDRH and seqCDRL.
    
Input: 
    seqCDRH, seqCDRL, they are the output of function SeqCDR
    
Output: 
    hcluster_CDRL, hcluster_CDRH. They are given in the form of
         [[0, 1, 0.5, 2],
          [2, 3, 0.5, 2],
          [4, 5, 1, 4]]
    each row is in the form of [idx1, idx2, dist, sample_count]
    In the above example, there are four samples with index 0, 1, 2, 3.
    [0, 1, 0.5, 2]: mean sample 0 and sample 1 are clustered together, them distance
    between sample 0 and sample 1 is 0.5. After clustering, there are 2 samples in this 
    cluster: sample 0 and sample 1.
    4 means the cluster formed in row number 4 - len(samples) = 0
    5 means the cluster formed in row number 5 - len(samples) = 1.
    Therefore, the last row means combine the clusters formed row 0 and row 1
    and the sample count is 4, which means all the samples(0, 1, 2, 3) are in this cluster.

'''

def Hcluster(seqCDRH, seqCDRL):

    from Bio import Align
    aligner = Align.PairwiseAligner()
    aligner.open_gap_score = 0 # it is unlikly that two sequence of the same source could be different by a gap
    aligner.extend_gap_score = 0
    aligner.match = 1
    aligner.mismatch = 0 
    aligner.mode = 'local'
    
    from scipy.cluster.hierarchy import linkage, optimal_leaf_ordering
    from scipy.spatial.distance import squareform
    import numpy as np
    
    CDR = [seqCDRH, seqCDRL]
    
    for idx in range(2):
        
        toydata = CDR[idx]
        distmatrix = np.zeros((len(toydata), len(toydata)))
        for i in range(len(toydata)):
            for j in range(len(toydata)):
                distmatrix[i,j] += aligner.score(toydata[i][1], toydata[j][1])
                distmatrix[i,j] += aligner.score(toydata[i][2], toydata[j][2])
                distmatrix[i,j] += aligner.score(toydata[i][3], toydata[j][3])
                l1 = len(toydata[i][1]) + len(toydata[i][2]) + len(toydata[i][3])
                l2 = len(toydata[j][1]) + len(toydata[j][2]) + len(toydata[j][3])
                l = min(l1, l2)
                distmatrix[i,j] =  1 - distmatrix[i,j]/l
                
        if idx == 0:
            dm_CDRH = squareform(distmatrix)
            hcluster_CDRH_unordered = linkage(dm_CDRH, method = 'complete')
            hcluster_CDRH = optimal_leaf_ordering(hcluster_CDRH_unordered, dm_CDRH )
        elif idx == 1:
            dm_CDRL = squareform(distmatrix)
            hcluster_CDRL_unordered = linkage(dm_CDRL, method = 'complete')
            hcluster_CDRL = optimal_leaf_ordering(hcluster_CDRL_unordered, dm_CDRL )
            
    return hcluster_CDRH, hcluster_CDRL

'''
Cut_distance_n_cluster:
    A function to calculate the number of clusters under different cut distances

Inputs:
    hcluster_CDRH, hcluster_CDRL: are the outputs of function Hcluster
    sd: the saving directory
Outputs:
    heights_nCDRH_nCDRL:
        a list in the form of:
            [[0.1, 0.2,...], [10, 5, ...], [8, 6, ...]]
        where [0.1, 0.2,...] gives the cut_distance, [10, 5, ...] gives the number of 
        clusters of CDRH if all the CDRHs are clustered togher if the distances are 
        no larger than the given distances. [8, 6, ...] are the cluster numbers of CDRL, 
        and interpreted similarly as CDRH.
        
        at cut distance 0.1, if all the samples are clustered together if the distance 
        is <= 0.1, then, CDRHs can be clustered into 10 clusters.
    
'''
                
def Cut_distance_n_cluster( hcluster_CDRH, hcluster_CDRL):
    
    from scipy.cluster.hierarchy import cut_tree
    
    n_clusters_CDRH = []
    n_clusters_CDRL = []
    heights = [x * 0.02 for x in range(51)]
    for h in heights:
        # Cut the tree at different heights
        cut_CDRH = cut_tree(hcluster_CDRH, height = h)
        cut_CDRL = cut_tree(hcluster_CDRL, height = h)
        # Calculate the number of clusters at different heights
        maxH = -1
        maxL = -1
        for i in range(len(cut_CDRH)):
            if cut_CDRH[i][0] >= maxH:
                maxH = cut_CDRH[i][0]
        n_clusters_CDRH.append(maxH)
        
        for i in range(len(cut_CDRL)):
            if cut_CDRL[i][0] >= maxL:
                maxL = cut_CDRL[i][0]
        n_clusters_CDRL.append(maxL)
    
    heights_nCDRH_nCDRL = [heights, n_clusters_CDRH, n_clusters_CDRL]
    
    return heights_nCDRH_nCDRL

'''
Draw_elbow:
    A function to draw the curves between the cut_distances and the number of clusters 
    of CDRH and CDRL. 
Input:
    heights_nCDRH_nCDRL: 
        the output of Cut_distance_n_cluster
    sd:
        the saving directory of the figures
Output:
    Save the figures to the working directory
'''
def Draw_elbow(sd, heights_nCDRH_nCDRL):
    os.chdir(sd)
    # Take out the values
    heights = heights_nCDRH_nCDRL[0]
    n_clusters_CDRH = heights_nCDRH_nCDRL[1]
    n_clusters_CDRL = heights_nCDRH_nCDRL[2]
    
    from matplotlib import pyplot as plt
    
    plt.title('CDRH_clusters Vs cut_distance')
    plt.xlabel('cut_distance')
    plt.ylabel('CDRH_clusters')
    plt.plot(heights, n_clusters_CDRH)
    plt.savefig('CDRH_clusters Vs cut_distance.png')
    plt.show()
    plt.close()
    
    plt.title('CDRL_clusters Vs cut_distance')
    plt.xlabel('cut_distance')
    plt.ylabel('CDRL_clusters')
    plt.plot( heights, n_clusters_CDRL)
    plt.savefig('CDRL_clusters Vs cut_distance.png')
    plt.show()
    plt.close()
################################################################################    
   
'''
Cluster_by_cut_dist:
    A function to find the clusters according to the giving cut_dist
Input:
    hcluster:
        one of the returned values of the Hcluster
    cut_dist:
        the cut distance
Output:
    a list of list, each list gives the indeces of the samples in this cluster.
'''              
def Cluster_by_cut_dist(hcluster, cut_dist):

    from scipy.cluster.hierarchy import cut_tree
    
    cut_CDR = cut_tree(hcluster, height = cut_dist)

            
    #find the corresponding ids
    n = 0
    for i in range(len(cut_CDR)):
        if cut_CDR[i][0] >= n:
            n = cut_CDR[i][0]

    clusters = []
    for i in range(n+1):
        subcluster = []
        for j in range(len(cut_CDR)):
            if cut_CDR[j][0] == i:
                subcluster.append(j)
        clusters.append(subcluster)
            
    return clusters
########################################################################
'''
Select_representatives:
    this function is to select the representive element, which has the
    largest contact number, in each cluster
Input:
    contact: 
        the output of AAC_2
    seqCDR:
        the same as above
    clusters:
        the output of Cluster_by_cut_dist
    
    ###the clusters and the seqCDR should be corresponding with each other, that
        is they are both about the CDRH or CDRL.###
Output:
    rpts:
        a list of the representative elements
            
'''

def Select_representatives(contact, seqCDR, clusters):
        # select the representatives 
    rpts = []
#    m_ctns = []
    for cluster in clusters:
        representative = -1
        max_contact = -1
        for i in cluster:
            contact_number = 0
            pdbid = seqCDR[i][0][:4]
            hl = seqCDR[i][0][4]
            HL = seqCDR[i][0][5]
            for fcdn in contact[pdbid]:
                if fcdn[0][0] == hl and fcdn[0][2] == HL:
                    contact_number += fcdn[3]
                    
            if contact_number >= max_contact:
                max_contact = contact_number
                representative = i
        rpts.append(seqCDR[representative][0])
#        m_ctns.append(max_contact)
            
    return rpts
##########################################################################
'''
Prepare_for_FrameConstraint:
    The output of this function is used as the beginning data of FrameConstraint module
Input:
    contact: the same as above, which is the one of the returned values of AAC_2
    CDRH_rpts: 
        the output of Select_representatives corresponding to CDRH
    CDRL_rpts:
        similar to CDRH_rpts
Output:
    ac_contact:
        it is the beginning data of FrameConstraint module.
'''                
def Prepare_for_FrameConstraint(contact, CDRH_rpts, CDRL_rpts):
    ac_contact = {}
    all_rep = copy.deepcopy(CDRH_rpts)
    all_rep.extend(CDRL_rpts)
    for rpt in all_rep:
        container = []
        pdbid = rpt[:4]
        hl = rpt[4]
        HL = rpt[5]
        for fcdn in contact[pdbid]:
            if fcdn[0][0] == hl and fcdn[0][2] == HL:
                container.append(fcdn)
        ac_contact[pdbid + HL + hl] = container
        
    return ac_contact
###########################################################################

def main(): 
    # wd is the working directory, where you can find the returned results from AAC_2
    # sd is the saving directory               
    wd= "/home/leo/Documents/Database/Pipeline_New/All with peptide 5+ resolution 4A"
    sd = '/home/leo/Documents/Database/Pipeline_New/Results'

    os.chdir(wd)
    # We have to process the training and the testing data separetely
    train_test = ['training', 'testing']
    for header in train_test:
        with open(header+'_sequence', 'r') as f:
            sequence = json.load(f)
        with open(header+'_matched_ids', 'r') as f:
            matched_ids = json.load(f)
        with open(header+'_contact', 'r') as f:
            contact = json.load(f)  
        # Do the alignment and selection   
        seqCDRH, seqCDRL = SeqCDR(sequence, matched_ids)
        
        hcluster_CDRH, hcluster_CDRL = Hcluster(seqCDRH, seqCDRL)
        
        clusters_CDRL = Cluster_by_cut_dist(hcluster_CDRL, 0.1)
        clusters_CDRH = Cluster_by_cut_dist(hcluster_CDRH, 0.1)
        
        CDRL_rpts = Select_representatives(contact, seqCDRL, clusters_CDRL)
        CDRH_rpts = Select_representatives(contact, seqCDRH, clusters_CDRH)
        
        ac_contact = Prepare_for_FrameConstraint(contact, CDRH_rpts, CDRL_rpts)
        
        heights_nCDRH_nCDRL = Cut_distance_n_cluster( hcluster_CDRH, hcluster_CDRL)
        
        # Save the results
        os.chdir(sd)
        Draw_elbow(sd, heights_nCDRH_nCDRL)
        with open(header+'_heights_nCDRH_nCDRL', 'w') as f:
            json.dump(heights_nCDRH_nCDRL, f)
        with open(header + 'ac_contact', 'w') as f:
            json.dump(ac_contact, f)
############################################################################
    
if __name__ == '__main__':
    main()


#working_directory = ["/home/leo/Documents/Database/Pipeline/All with peptide 5+ resolution 4A",\
#                     '/home/leo/Documents/Database/Pipeline/Complexes with Affinity']
#os.chdir(working_directory[1])
#
#with open('sequence', 'r') as f:
#    sequence = json.load(f)
#with open('good_matched_ids', 'r') as f:
#    matched_ids = json.load(f)
#with open('contact', 'r') as f:
#    contact = json.load(f)        
#
#seqCDRH, seqCDRL = SeqCDR(sequence, matched_ids)
#
#hcluster_CDRH, hcluster_CDRL = Hcluster(seqCDRH, seqCDRL)
#
#clusters_CDRL = Cluster_by_cut_dist(hcluster_CDRL, 0.1)
#clusters_CDRH = Cluster_by_cut_dist(hcluster_CDRH, 0.1)
#CDRL_rpts = Select_representatives(contact, seqCDRL, clusters_CDRL)
#CDRH_rpts = Select_representatives(contact, seqCDRH, clusters_CDRH)
#ac_contact = Prepare_for_FrameConstraint(contact, CDRH_rpts, CDRL_rpts)
#heights_nCDRH_nCDRL = Cut_distance_n_cluster( hcluster_CDRH, hcluster_CDRL)
#sd = '/home/leo/Documents/Database/Pipeline_New/Results'
#Draw_elbow(sd, heights_nCDRH_nCDRL)

