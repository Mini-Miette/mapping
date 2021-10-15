#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 09:59:33 2021

@author: lxenard
"""

from collections import Counter
import pandas as pd
import time


def get_bwt(sequence):
    """ Calcule la transformée de Burrows Wheeler.

    Args:
        sequence (str): La séquence à transformer.

    Returns:
        str: Transformée de Burrows Wheeler.

    """
    sequence += '$'
    permutations = [sequence]
    for i in range(len(sequence)-1):
        permutations.append(permutations[i][-1] + permutations[i][:-1])
    permutations.sort()
    return ''.join([permut[-1] for permut in permutations])


def get_p_n(bwt):
    """


    Args:
        bwt (str): Transformée de Burrows Wheeler.

    Returns:
        p (list(int)): DESCRIPTION.
        n (dict(str : int)): Position de la première occurrence des lettres.

    """

    alphabet = list(set(bwt))
    alphabet.sort()

    counts = []
    for letter in alphabet:
        counts.append(bwt.count(letter))

    p = [bwt[:i].count(letter) for i, letter in enumerate(bwt)]

    n = {'$': 0}
    for i in range(1, len(alphabet)):
        previous_pos = n[alphabet[i-1]]
        n[alphabet[i]] = counts[i-1] + previous_pos

    return p, n


def get_reverse_bwt(bwt):
    """ Calcule la transformée inverse de Burrows Wheeler.

    Args:
        bwt (str): Transformée de Burrows Wheeler.

    Returns:
        str: Transformée inverse de Burrows Wheeler.

    """
    p, n = get_p_n(bwt)

    # Initialisation de la séquence avec des 0 et on va la remplir
    # petit à petit.
    ori_seq = [0 for element in bwt]
    # Ajoute des positions connues.
    dollar_index = bwt.find('$')
    ori_seq[-1] = '$'
    ori_seq[-2] = bwt[dollar_index - 1]

    # Récupération de la lettre et du p courants.
    current_p = p[dollar_index - 1]
    current_letter = ori_seq[-2]

    # Itération sur les positions vides de la séquence.
    for i_seq in range(len(bwt) - 3, -1, -1):
        new_letter_index = n[current_letter] + current_p
        ori_seq[i_seq] = bwt[new_letter_index]

        # Mise à jour de la lettre et du p courants.
        current_p = p[new_letter_index]
        current_letter = bwt[new_letter_index]

    return ''.join(ori_seq[:-1])


def get_reverse_bwt_position(bwt):
    """
    Va servir à vérifier que la position du read est correcte.

    Args:
        bwt (TYPE): DESCRIPTION.

    Returns:
        pos (TYPE): DESCRIPTION.

    """
    p, n = get_p_n(bwt)
    i = len(bwt) - 2
    pos = [-1] * len(bwt)
    j = 0
    while bwt[j] != '$':
        pos[j] = i
        j = n[bwt[j]] + p[j]
        i -= 1

    return pos


def build_FM_index(bwt):

    alphabet = list(set(bwt))
    alphabet.sort()

    df = pd.DataFrame([{letter: 0 for letter
                        in alphabet}] * (len(bwt) + 1))

    for i in range(1, len(bwt)+1):
        letter = bwt[i-1]
        df[letter][i:] += 1

    return df


def match(read, bwt, N, FM):
    """
    return d et f

    Args:
        read (TYPE): DESCRIPTION.
        bwt (TYPE): DESCRIPTION.
        N (TYPE): DESCRIPTION.
        FM (TYPE): DESCRIPTION.

    Returns:
        None.

    """
    alphabet = list(set(bwt))
    alphabet.sort()

    letter = read[-1]
    begin = N[letter]
    end = N[letter] + FM[letter][len(bwt)] - 1

    for letter in read[-2::-1]:
        begin = N[letter] + FM[letter][begin]
        end = N[letter] + FM[letter][end + 1] - 1

        if end < begin:
            begin = -1
            end = -1
            break

    return (begin, end)


def read_fasta(fasta_file):
    """


    Args:
        fasta_file (TYPE): DESCRIPTION.

    Yields:
        sequence (str): DESCRIPTION.
        reference (str): DESCRIPTION.

    """

    with open(fasta_file, 'r') as my_file:

        activeone = False
        sequence = ''
        reference = ''
        for line in my_file:
            if str(line).startswith(">"):
                if activeone:
                    yield (sequence, reference)
                else:
                    activeone = True
                reference = line.strip()[1:]
                line = next(my_file, None)
                sequence = str(line).strip()
            else:
                sequence += str(line).strip()
        yield (sequence, reference)


def get_pos(begin, end, pos_list):
    if begin == -1 and end == -1:
        position = [-1]
    else:
        position = []
        for i in range(begin, end + 1):
            position.append(pos_list[i] + 1)
    return position


if __name__ == '__main__':

    start_time = time.time()

# =============================================================================
#     fasta_file = '/home/sdv/m2bi/lxenard/Documents/Omiques/test_ref.fasta'
#     read_file = '/home/sdv/m2bi/lxenard/Documents/Omiques/test_read.fasta'
#     output_file = '/home/sdv/m2bi/lxenard/Documents/Omiques/test_output.bi'
# =============================================================================

# =============================================================================
#     fasta_file = '/home/sdv/m2bi/lxenard/Documents/Omiques/masterBI/NC_045512-N.fna'
#     read_file = '/home/sdv/m2bi/lxenard/Documents/Omiques/masterBI/reads_1000_10_patient_6.fna'
#     output_file = '/home/sdv/m2bi/lxenard/Documents/Omiques/match_1000_10_patient_6.bi'
# =============================================================================

# =============================================================================
#     fasta_file = '/home/sdv/m2bi/lxenard/Documents/Omiques/masterBI/NC_045512-N.fna'
#     read_file = '/home/sdv/m2bi/lxenard/Documents/Omiques/masterBI/reads_10000_30_patient_6.fna'
#     output_file = '/home/sdv/m2bi/lxenard/Documents/Omiques/match_10000_30_patient_6.bi'
# =============================================================================

    fasta_file = '/home/sdv/m2bi/lxenard/Documents/Omiques/masterBI/NC_045512-N.fna'
    read_file = '/home/sdv/m2bi/lxenard/Documents/Omiques/masterBI/reads_10000_100_patient_6.fna'
    output_file = '/home/sdv/m2bi/lxenard/Documents/Omiques/match_10000_100_patient_6.bi'

    # Lecture de la séquence de référence.
    ref_reader = read_fasta(fasta_file)
    ref_seq, ref_name = next(ref_reader)
    bwt = get_bwt(ref_seq)
    P, N = get_p_n(bwt)
    FM = build_FM_index(bwt)

    reader = read_fasta(read_file)
    for (read_seq, read_name) in reader:
        begin, end = match(read_seq, bwt, N, FM)
        pos_list = get_reverse_bwt_position(bwt)
        positions = get_pos(begin, end, pos_list)

        with open(output_file, 'a') as output:
            for pos in positions:
                output.write(f'{ref_name}\t{read_name}\t{pos}\n')

    total_time = time.time() - start_time
    print('\nDONE in {:.0f} min {:.2f} s.'.format(total_time // 60,
                                                  total_time % 60))

    # Récupérer :
    #   - #
    #   - 0 match
    #   - 1 match
    #   - 2+ matches

    with open(output_file, 'r') as file:
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


# =============================================================================
# p, n = get_p_n(bwt)
# i = len(bwt) - 1
# orig = ""
# j = 0
# while i >= 1:
#     print(orig)
#     orig = bwt[j] + orig
#     j = n[bwt[j]] + p[j]
#     i -= 1
# =============================================================================
