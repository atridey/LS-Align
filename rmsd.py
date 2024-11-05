from Bio.PDB import PDBParser
from math import pow, sqrt
#7.550 w/ mobile pdb, target gro
#Training 4 mobile 5 target, 0.008


#Finds RMSD given PDB files
#Input: file1 (String, filename 1), file2 (String, filename 2)
#Output: Float, root mean square deviation
def rmsd(file1, file2):
    #Finding coordinates
    temp_struct = PDBParser(QUIET = True).get_structure(file1, file1)
    coords1 = {}
    for model in temp_struct:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    coords1.update({atom.get_name(): atom.get_coord()})
    temp_struct = PDBParser(QUIET = True).get_structure(file2, file2)
    coords2 = {}
    for model in temp_struct:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    coords2.update({atom.get_name(): atom.get_coord()})
    
    #Calculation
    sum = 0
    n = 0
    for i in coords1:
        try:
            sum += pow(coords2[i][0]-coords1[i][0], 2) + pow(coords2[i][1]-coords1[i][1], 2) + pow(coords2[i][2]-coords1[i][2], 2) #Square root and square cancel out, so nothing else is done
            n += 1
        except:
            pass
    return sqrt((1/n) * sum)


file1 = input('Filename 1: ')
file2 = input('Filename 2: ')
print(rmsd(file1, file2))
