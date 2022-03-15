from rdkit import Chem
import os
from pathlib import Path

def findcharge(sdffile) -> list:
    with open(sdffile, "r") as f:
        for line in f.readlines():
            if line.startswith("M  CHG"):
                line = line.strip()
                line = line.split()
                atom_num = int(line[2])
                idx = 3
                charged_idx = []
                for i in range(atom_num):
                    charged_idx.append(int(line[idx]))
                    idx += 2
        # charged_atom is now a list of all charged_atoms in mol
    return charged_idx

def getsdfchg(sdffile):

    sdflines = Path(sdffile).read_text().split("\n")
    ifchg = False
    for line in sdflines:
        if line.startswith("M  CHG"):
            ifchg = True
            line = line.strip()
            line = line.split()
            atom_num = int(line[2])
            counter = 4
            charge = []
            for i in range(atom_num):
                charge.append(int(line[counter]))
                counter +=2
            assert len(charge) == atom_num, 'something wrong'
            total_charge = sum(charge)

    chg = 0 if not ifchg else total_charge

    return chg

def findneighbour(mol, idxr) -> list:
    atom = mol.GetAtomWithIdx(idxr)  # Atoms: c(harged)atom
    neighbours = atom.GetNeighbors()  # list of Atoms: neighbour atoms of catom
    neigh_idxr = [
        x.GetIdx() for x in neighbours
    ]  # list of int: idxes of neighbours, include hydrogens
    # neigh_idxr = [x.GetIdx() for x in neighbours if x.GetAtomicNum != 1] #list of int: idxes of neighbours, without hydrogens

    return neigh_idxr


def neutralize_atoms(mol):
    # Neutralize molecules for each charged atom
    pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
    at_matches = mol.GetSubstructMatches(pattern)
    at_matches_list = [y[0] for y in at_matches]
    if len(at_matches_list) > 0:
        for at_idx in at_matches_list:
            atom = mol.GetAtomWithIdx(at_idx)
            chg = atom.GetFormalCharge()
            hcount = atom.GetTotalNumHs()
            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(hcount - chg)
            atom.UpdatePropertyCache()
    return mol


def GetRingSystems(mol, includeSpiro=False):  # do not count Spiro into account
    # High memory usage
    ri = mol.GetRingInfo()
    systems = []
    for ring in ri.AtomRings():
        ringAts = set(ring)
        nSystems = []
        for system in systems:
            nInCommon = len(ringAts.intersection(system))
            if nInCommon and (includeSpiro or nInCommon > 1):
                ringAts = ringAts.union(system)
            else:
                nSystems.append(system)
        nSystems.append(ringAts)
        systems = nSystems
    return systems

def Get_sorted_heavy(mol, idxlist):
    atoms = [ mol.GetAtomWithIdx(idx) for idx in idxlist ]
    atom_nums = [atom.GetAtomicNum() for atom in atoms]

    for element in [17,16,15,9,8,7,6]: #Cl,S,P,F,O,N,C
        if element in atom_nums: 
            return idxlist[atom_nums.index(element)]
    
    return None # if no selected atoms

# Read mol2 molecule once per time
def next_mol2_lines(infile):
    """Method to return one mol2 block once."""
    lines = list()

    for line in open(infile):
        if "@<TRIPOS>MOLECULE" in line:
            if len(lines) == 0:
                lines.append(line)
            else: # in case there are multiple mol2blocks in infile
                yield lines
                lines = list()
                lines.append(line)
        else:
            lines.append(line)

    yield lines
