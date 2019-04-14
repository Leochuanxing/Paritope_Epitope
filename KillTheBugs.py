import os
import json

os.chdir('/home/leo/Documents/Database/Pipeline_New/Latest')
summary = open('summary.tsv', 'r')# case sensitive
file = summary.readlines()

good_ids = []
for line in file:
    line_list = line.split('\t')
    if line_list[0] == '6mid':
        print(line_list)
    if len(line_list) > 5:
        antigen_type = line_list[5].split('|')
        for st in antigen_type:
            st.strip()
        if 'peptide' in antigen_type or 'protein' in antigen_type:
            good_ids.append(line_list[0])
good_ids

len(good_ids)

with open('matched_ids', 'r') as f:
    matched_ids = json.load(f)
with open('combined_ids', 'r') as f:
    combined_ids = json.load(f)
with open('sequence', 'r') as f:
    sequence = json.load(f)
    


