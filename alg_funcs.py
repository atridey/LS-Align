from openbabel import pybel
from openbabel.openbabel import *
from math import pow
import numpy as np
np.set_printoptions(suppress = True)
from scipy.spatial.transform import Rotation


#Input: query (pybel molecule), template (pybel molecule)
#Output: scores (tuple with mass, radius, bond), alignments (tuple with mass, radius, bond)
def greedy_search(query, template):
    #Mass
    score_mass = np.empty([len(query.atoms), len(template.atoms)])
    for query_index in range(len(query.atoms)): 
        for template_index in range(len(template.atoms)):
            score_mass[query_index, template_index] = (1/(1+ pow( (query.atoms[query_index].atomicmass - template.atoms[template_index].atomicmass), 2) / 100))
    #Radius
    score_radius = np.empty([len(query.atoms), len(template.atoms)])
    for query_index in range(len(query.atoms)): 
        for template_index in range(len(template.atoms)):
            score_radius[query_index, template_index] = (1/(1+ pow( (GetVdwRad(query.atoms[query_index].atomicnum) - template.atoms[template_index].atomicnum), 2) / 1.44))
    #Bond
    bond_weights = {'1': 1, '2': 1.2, '3': 1.4, 'am': 1.8, 'ar': 2, 'du': 1.6}
    mol2 = query.write('mol2').splitlines()
    bonds_line = None
    for line in range(len(mol2)):
        if (mol2[line] == "@<TRIPOS>BOND"):
            bonds_line = line+1
            break
    mol2 = mol2[bonds_line:]
    query_bonds = []
    for bond in mol2:
        query_bonds.append(bond.split())
    mol2 = template.write('mol2').splitlines()
    bonds_line = None
    for line in range(len(mol2)):
        if (mol2[line] == "@<TRIPOS>BOND"):
            bonds_line = line+1
            break
    mol2 = mol2[bonds_line:]
    template_bonds = []
    for bond in mol2:
        template_bonds.append(bond.split())
    score_bond = np.empty([len(query.atoms), len(template.atoms)])
    for query_index in range(len(query.atoms)):
        for template_index in range(len(template.atoms)):
            qatom_bonds = []
            for i in range(len(query_bonds)):
                current_bond = query_bonds[i]
                if (query_index+1 == int(current_bond[1])) or (query_index+1 == int(current_bond[2])):
                    qatom_bonds.append(current_bond[3])
            tatom_bonds = []
            for i in range(len(template_bonds)):
                current_bond = template_bonds[i]
                if (template_index+1 == int(current_bond[1])) or (template_index+1 == int(current_bond[2])):
                    tatom_bonds.append(current_bond[3])
            overlap = tuple(set(qatom_bonds).intersection(set(tatom_bonds)))
            qatom_weighted = 0
            for bond in qatom_bonds:
                qatom_weighted += bond_weights[bond]
            tatom_weighted = 0
            for bond in tatom_bonds:
                tatom_weighted += bond_weights[bond]
            overlap_weighted = 0
            for bond in overlap:
                overlap_weighted += bond_weights[bond]
            score_bond[query_index, template_index] = (overlap_weighted / (qatom_weighted + tatom_weighted - overlap_weighted) )
    scores = (score_mass, score_radius, score_bond)

    #Greedy search - Alignment
    #Mass
    score_copy = score_mass.copy()
    align_mass = []
    while (score_copy.max() != 0):
        max_coords = np.unravel_index(np.argmax(score_copy), score_copy.shape)
        align_mass.append( (max_coords[0], max_coords[1]) )
        score_copy[max_coords[0], :] = 0
        score_copy[:, max_coords[1]] = 0
    #Radius
    score_copy = score_radius.copy()
    align_radius = []
    while (score_copy.max() != 0):
        max_coords = np.unravel_index(np.argmax(score_copy), score_copy.shape)
        align_radius.append( (max_coords[0], max_coords[1]) )
        score_copy[max_coords[0], :] = 0
        score_copy[:, max_coords[1]] = 0
    #Bond
    score_copy = score_bond.copy()
    align_bond = []
    while (score_copy.max() != 0):
        max_coords = np.unravel_index(np.argmax(score_copy), score_copy.shape)
        align_bond.append( (max_coords[0], max_coords[1]) )
        score_copy[max_coords[0], :] = 0
        score_copy[:, max_coords[1]] = 0
    alignments = (sorted(align_mass), sorted(align_radius), sorted(align_bond))

    return scores, alignments

#Input: score (ndarray), alignment (list of tuples)
#Output: alignment (list of tuples)
def enhanced_greedy_search(score, alignment):
    while True:
        changed = False
        delta_score = []
        pairs = []

        for pair in alignment:
            i = pair[0]
            j = pair[1]
            for column_index in range(score.shape[1]):
                if column_index == j:
                    continue
                else:
                    l = column_index

                for row_index in range(score.shape[0]):
                    if (row_index, l) in alignment:
                        k = row_index
                        delta_score.append(  (score[i,l]-score[i,j]) + (score[k,j]-score[k,l])  )
                        pairs.append((i,j,k,l))
                        break
                    
        s_max = max(delta_score)
        if s_max > 0:
            changed = True
            print(s_max)
            pairs = pairs[delta_score.index(s_max)]
            index_1 = alignment.index((pairs[0], pairs[1])) #Find index old i,j
            index_2 = alignment.index((pairs[2], pairs[3])) #Old k,l
            alignment[index_1] = (pairs[0], pairs[3]) #(i,j)-->(i,l)
            alignment[index_2] = (pairs[2], pairs[1]) #(k,l)-->(k,j)

        if changed == False:
            break
    return sorted(alignment)


#Input: query (pybel molecule), template (pybel molecule), alignment (list of alignment pairs)
#Output: Aligned query vectors (2D ndarray)
def kabsch_alignment(query, template, alignment):
    #Apply initial alignment
    query_vectors = np.empty([len(query.atoms), 3])
    for pair in alignment:
        query_vectors[pair[1]] = query.atoms[pair[0]].coords
    template_vectors = []
    for atom in template.atoms:
        template_vectors.append(atom.coords)
    rotation = Rotation.align_vectors(template_vectors, query_vectors)[0] #Transform b onto a, then get rid of RSSD
    return rotation.apply(query_vectors)

#Input: template (pybel molecule), aligned_vectors (2D ndarray)
#Output: scores (2D ndarray)
def ls_score(template, aligned_vectors):
    #d0 = 2.0
    scores = np.empty([len(aligned_vectors), len(template.atoms)])
    for i in range(len(aligned_vectors)):
        for j in range(len(template.atoms)):
            i_coords = aligned_vectors[i]
            j_coords = np.array(template.atoms[j].coords)
            scores[i,j] = (1/(1+ pow(np.linalg.norm(i_coords-j_coords), 2) / 4) )
    return scores

#Input: template (pybel molecule), aligned_vectors (2D ndarray), alignment (list of tuples)
#Output: score (int)
def mol_ls_score(aligned_vectors, template, alignment):
    #d0 = 2.0
    sum = 0
    for pair in alignment:
        i_coords = aligned_vectors[pair[0]]
        j_coords = template.atoms[pair[1]].coords
        distSq = pow( np.linalg.norm(i_coords-j_coords), 2 )
        sum += (1/ (1 + distSq / 4))
    return((1/len(template.atoms)) * sum)

def ls_trim(ls_matrix, lower_limit, upper_limit):
    for row in range(0, lower_limit):
        ls_matrix[row] = np.zeros(ls_matrix.shape[0])
    for row in range(lower_limit, upper_limit):
        ls_matrix[row][:lower_limit] = np.zeros(lower_limit)
        ls_matrix[row][upper_limit:] = np.zeros(ls_matrix.shape[0] - upper_limit)
    for row in range(upper_limit, ls_matrix.shape[1]):
        ls_matrix[row] = np.zeros(ls_matrix.shape[0])
    return ls_matrix

#Input: query (pybel molecule), template (pybel molecule), initial_alignment (list of tuples) , lower_limit (int), upper_limit (int)
#Output: alignment (list of tuples)
def single_iteration(query, template, initial_alignment, lower_limit, upper_limit):
    initial_alignment = initial_alignment[lower_limit:upper_limit]
    #Initial EGS (as if it's a do while loop)
    aligned_vectors = kabsch_alignment(query, template, initial_alignment)
    ls_matrix = ls_score(template, aligned_vectors)
    #ls_matrix = ls_trim(ls_matrix, lower_limit, upper_limit)
    old_egs_alignment = enhanced_greedy_search(ls_matrix, initial_alignment) #Very slow
    while True:
        aligned_vectors = kabsch_alignment(query, template, old_egs_alignment)
        ls_matrix = ls_score(template, aligned_vectors)
        #ls_matrix = ls_trim(ls_matrix, lower_limit, upper_limit)#np.array([ls_matrix[0, lower_limit:upper_limit], ls_matrix[1, lower_limit:upper_limit]])
        new_egs_alignment = enhanced_greedy_search(ls_matrix, old_egs_alignment) #Very slow
        if old_egs_alignment == new_egs_alignment:
            return new_egs_alignment

#Input: query (pybel molecule), template (pybel molecule), alignment (list of tuples), 
#Output: alignment (list of tuples)
def heuristic_iteration(query, template, alignment):
    alignment = single_iteration(query, template, alignment, 0, len(query.atoms))
    length = round(len(query.atoms) / 2)
    while (length >= 4):
        lower_limit = 0
        upper_limit = length
        while(upper_limit <= len(query.atoms)):
            print(length, lower_limit)
            align_frag = single_iteration(query, template, alignment, lower_limit, upper_limit)
            for i in range(lower_limit, upper_limit):
                alignment[i] = align_frag[i - lower_limit] 
            lower_limit += 1
            upper_limit += 1
        length = round(length / 2)
    return alignment

#Input: query (pybel molecule), template (pybel molecule), alignment (list of tuples), 
#Output: alignment (list of tuples)
def not_heuristic_iteration(query, template, alignment):
    alignment = single_iteration(query, template, alignment, 0, len(query.atoms))
    return alignment
