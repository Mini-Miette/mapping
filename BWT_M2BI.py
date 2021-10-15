#!/bin/python3

import sys

def gener_BWT_transform(text):
    """Burrows-Wheeler transformation
    """
    bwt_text = []
    for letter in text :
        text = text[-1:] + text[:-1]
        bwt_text.append(text)
    return "".join([i[-1] for i in sorted(bwt_text)])
 
 
def gener_P_N(transform):
    """N & P matrices for Burrows-Wheeler transformation
    """
    P = []
    dico = {"$" : 0, "A" : 0, "C" : 0, "G" : 0, "T" : 0}
    for lettre in transform :
        P.append(dico[lettre])
        dico[lettre] += 1
    N = {}
    count = 0
    for k in dico.keys() :
        N[k] = count
        count += dico[k]
    return P, N
 
def gener_retro_bwt(transform, P, N) :
    i = len(transform) - 1
    origin = "$"
    j = 0
    while i >= 1 :
        origin = transform[j] + origin
        j = N[transform[j]] + P[j]
        i -= 1
    return origin
 
def pos_BWT(transform, P, N):
    i = len(transform) - 2
    pos = [-1] * len(transform)
    j = 0
    while transform[j] != "$" :
        pos[j] = i
        j = N[transform[j]] + P[j]
        i -= 1
    return pos
 
def FM_index(transform):
    dict_list = []
    dicto = {"$" : 0, "A" : 0, "C" : 0, "G" : 0, "T" : 0}
    dict_list.append(dicto.copy())
    for letter in transform :
        dicto[letter] += 1
        dict_list.append(dicto.copy())
    return dict_list
 
def match(read, N, FM):
    nt = read[-1]
    d = N[nt] #debut
    f = N[nt] + FM[-1][nt] -1
    # print(d,f)
    # print(read[-2::-1]+'$')
    # print("N="+str(N))
    # print(FM)
    for nt in read[-2::-1]:
        #print(nt)
        d = N[nt] + FM[d][nt]
        f = N[nt] + FM[f+1][nt] -1
        #print("d-f",d,f)
        if f < d :
            #print('pouet')
            return d, f
    return d, f

def find_pos(start, end,  pos_list):
    position = []
    #print(start)
    #print(end)
    if start < end:
      return [-1]
    for i in range(start, end+1):
        position.append(pos_list[i]+ 1) 
    return position

all_args = sys.argv
#print(all_args)

file_ref = open(sys.argv[1], mode = 'r')
ref = file_ref.read()
ref = ref.split('\n')
ref_name = str(ref[0])
ref = ''.join(ref[1:])
file_ref.close()
ref = ref + '$'
bwt_ref = gener_BWT_transform(ref)
P_N = gener_P_N(bwt_ref)
fm_index = FM_index(bwt_ref)
pos = pos_BWT(transform = bwt_ref, P = P_N[0] ,N = P_N[1])
#print(ref_name)
#print(bwt_ref)


file_read = open(sys.argv[2], mode = 'r')

for line in file_read:
  line = line.rstrip()
  #print(line.rstrip())
  if (line.startswith('>')):
    curr_read = line
    continue
  curr_match = match(read = line, N = P_N[1], FM = fm_index)
  curr_pos = find_pos(start = curr_match[1], end = curr_match[0], pos_list = pos)
  
  for posi in curr_pos:
    output = f'{ref_name}\t{curr_read}\t{str(posi)}'
    print(output)