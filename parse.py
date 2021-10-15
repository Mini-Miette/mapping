#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 15 16:07:46 2021

@author: lxenard
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 09:59:33 2021

@author: lxenard
"""

# Récupérer :
#   - #
#   - 0 match
#   - 1 match
#   - 2+ matches

# Faire représentation graphique

#bi_file = '/home/sdv/m2bi/lxenard/Documents/Omiques/test_output.bi'
#bi_file = '/home/sdv/m2bi/lxenard/Documents/Omiques/match_1000_10_patient_6.bi'
from collections import Counter
bi_file = '/home/sdv/m2bi/lxenard/Documents/Omiques/match_10000_30_patient_6.bi'
#bi_file = '/home/sdv/m2bi/lxenard/Documents/Omiques/match_10000_100_patient_6.bi'

with open(bi_file, 'r') as file:
    match_count = Counter()
    no_match = []
    for line in file:
        ref, read, pos = line.split('\t')
        if pos.strip() == '-1':
            no_match.append(read)
        else:
            match_count[read] += 1

# print(match_count)

nb_no_match = len(no_match)
nb_1_match = Counter(match_count.values())[1]
nb_2_match = Counter(match_count.values())[2]

print(f'0 match: {nb_no_match}')
print(f'1 match: {nb_1_match}')
print(f'2+ match: {nb_2_match}')
