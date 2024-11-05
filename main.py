from openbabel import pybel
from openbabel.openbabel import *
from math import pow
import numpy as np
np.set_printoptions(suppress = True)
obErrorLog.SetOutputLevel(0)

from alg_funcs import *

#Reading
filename = input('Enter query file name: ')
filetype = filename.split('.')[1]
query = next(pybel.readfile(filetype, filename))

filename = input('Enter template file name: ')
filetype = filename.split('.')[1]
template = next(pybel.readfile(filetype, filename))

#GS
(score_mass, score_radius, score_bond), (align_mass, align_radius, align_bond) = greedy_search(query, template)
#EGS
align_mass = enhanced_greedy_search(score_mass, align_mass)
align_radius = enhanced_greedy_search(score_radius, align_radius)
align_bond = enhanced_greedy_search(score_bond, align_bond)
#Heuristic iteration
align_mass = heuristic_iteration(align_mass)
align_radius = heuristic_iteration(align_radius)
align_bond = heuristic_iteration(align_bond)

#Score and pick max
ls_mass = mol_ls_score(kabsch_alignment(query, template, align_mass), template, align_mass)
ls_radius = mol_ls_score(kabsch_alignment(query, template, align_radius), template, align_radius)
ls_bond = mol_ls_score(kabsch_alignment(query, template, align_bond), template, align_bond)

ls_max = max(ls_mass, ls_radius, ls_bond)

align_max = None
if ls_max == ls_mass:
    align_max = align_mass
elif ls_max == ls_radius:
    align_max = align_radius
else: #assuming equivalence to elif ls_max == ls_bond
    align_max = align_bond

for i, j in align_max:
    template.OBMol.AddHydrogens(False)
    query.OBMol.AddHydrogens(False)
    template.OBMol.SetAtomPDBIndex(i, 1)
    query.OBMol.SetAtomPDBIndex(j, 1)
        
    aligned_query = template.Align(query)

aligned_query.write('aligned_query.pdb')
print('Written to \"aligned_query.pdb\".')
