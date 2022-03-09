import os
from pickle import TRUE
import sys
import copy
from functools import partial
from pathlib import Path


from rdkit import Chem
from rdkit.Chem.Lipinski import RotatableBondSmarts
from rdkit.Chem.rdDistGeom import EmbedMolecule, EmbedMultipleConfs, ETKDGv3
from rdkit.Chem.rdForceFieldHelpers import MMFFOptimizeMolecule,MMFFOptimizeMoleculeConfs

sys.path.append("/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/utils")
from utils import GetRingSystems, findneighbour, Get_sorted_heavy
from testFrags import RDK_NCOPS_group_SMARTS_NOS_simplified

def checkbonded(mol, torsion_quartet):
    set_tor_quart = set([int(idxr) for idxr in torsion_quartet])
    for idxr in set_tor_quart:
        neigh_idxr2 = findneighbour(mol, int(idxr))
        if len(set(neigh_idxr2).intersection(set_tor_quart)) == 2:
            idxr2 = idxr
            break

    for idxr in set_tor_quart:
        neigh_idxr = findneighbour(mol, int(idxr))
        if len(set(neigh_idxr).intersection(set_tor_quart)) == 1:
            if idxr2 in neigh_idxr:
                idxr1 = idxr
            else:
                idxr4 = idxr
        else:
            if idxr != idxr2:
                idxr3 = idxr

    return idxr1, idxr2, idxr3, idxr4
    

def GetTorsion0(mol) -> list:
    TorsionPairs = []
    all_rot = mol.GetSubstructMatches(RotatableBondSmarts)
    for rotpair in all_rot:
        # idx2, idx3 = rotpair[2], rotpair[3]
        idx2, idx3 = rotpair[0], rotpair[1]

        atom2, atom3 = mol.GetAtomWithIdx(idx2), mol.GetAtomWithIdx(idx3)
        if atom2.GetAtomicNum() != 1 and atom3.GetAtomicNum() != 1:
            TorsionPairs.append([idx2, idx3])
    return TorsionPairs


def GetQuartetAtoms1(mol, TorsionPair):
    idxr2, idxr3 = TorsionPair
    nb2 = findneighbour(mol, idxr2)
    nb2.remove(idxr3)
    nb3 = findneighbour(mol, idxr3)
    nb3.remove(idxr2)

    nb2_noH = [i for i in nb2 if mol.GetAtomWithIdx(i).GetAtomicNum() != 1]
    nb3_noH = [i for i in nb3 if mol.GetAtomWithIdx(i).GetAtomicNum() != 1]
    
    # if nb2_noH and nb3_noH:
        # idxr1, idxr4 = nb2_noH[0], nb3_noH[0]
    idxr1, idxr4 = Get_sorted_heavy(mol, nb2), Get_sorted_heavy(mol, nb3) # only accept Cl,S,P,F,O,N,C, hetero prioritize， None otherwise
    if idxr1 != None and idxr4 != None:
        TorsionQuartet = [idxr1, idxr2, idxr3, idxr4]
        return TorsionQuartet
    else:
        pass

def ExpandQuartet2(mol, TorsionQuartet):
    idxr1, idxr2, idxr3, idxr4 = TorsionQuartet
    func = partial(findneighbour, mol)
    nb1, nb2, nb3, nb4 = map(func, [idxr1, idxr2, idxr3, idxr4])
    ExpandedTorsion = list(set(nb1 + nb2 + nb3 + nb4))

    return list(ExpandedTorsion)


def KeepRing3(mol, ExpandedTorsion, rings):
    for ring in rings:
        for idx in ExpandedTorsion:
            if idx in ring:
                ExpandedTorsion.extend(list(ring))
                break
    ExpandedTorsion = set(ExpandedTorsion)

    return list(ExpandedTorsion)


def KeepOrtho4(mol, ExpandedTorsion, TorsionQuartet, rings):
    rings = GetRingSystems(mol, includeSpiro=False)
    # rdkit.Chem.rdmolops.GetShortestPath((Mol)arg1, (int)arg2, (int)arg3) → tuple
    for ring in rings:
        for ring_idx in ring:
            if ring_idx in ExpandedTorsion:
                neighs = findneighbour(mol, ring_idx)
                label = False
                for neigh1 in neighs:
                    if neigh1 not in ring and neigh1 in TorsionQuartet:
                        label = True
                if label:
                    for neigh1 in neighs:
                        if neigh1 in ring: # in case of aliphatic ring, neigh1 must in the same ring
                            neigh2s = findneighbour(mol, neigh1)
                            ExpandedTorsion.extend(neigh2s)

    ExpandedTorsion = set(ExpandedTorsion)
    return list(ExpandedTorsion)


def KeepFuncGroup5(mol, ExpandedTorsion, FuncGroupIdxs):
    for funcgroup in FuncGroupIdxs: 
        for idx in ExpandedTorsion:
            # if mol.GetAtomWithIdx(idx).GetAtomicNum() not in [1,6]:
            if idx in funcgroup:
                ExpandedTorsion.extend(list(funcgroup))
                break
    
    ExpandedTorsion = set(ExpandedTorsion)

    ############## Not used now #############
    """
    # # adapted from lfg.py in rdkit contributions

    # atoms connected by non-aromatic double or triple bond to any heteroatom
    # c=O should not match (see fig1, box 15).  I think using A instead of * should sort that out?
    PATT_DOUBLE_TRIPLE = Chem.MolFromSmarts('A=,#[!#6]')
    # atoms in non aromatic carbon-carbon double or triple bonds
    PATT_CC_DOUBLE_TRIPLE = Chem.MolFromSmarts('C=,#C')
    # acetal carbons, i.e. sp3 carbons connected to tow or more oxygens, nitrogens or sulfurs; these O, N or S atoms must have only single bonds
    PATT_ACETAL = Chem.MolFromSmarts('[CX4](-[O,N,S])-[O,N,S]')
    # all atoms in oxirane, aziridine and thiirane rings
    PATT_OXIRANE_ETC = Chem.MolFromSmarts('[O,N,S]1CC1')

    PATT_TUPLE = (PATT_DOUBLE_TRIPLE, PATT_CC_DOUBLE_TRIPLE, PATT_ACETAL, PATT_OXIRANE_ETC)

    marked = []
    for patt in PATT_TUPLE:
        for path in mol.GetSubstructMatches(patt):
            for atomindex in path:
                marked.add(atomindex)
    """
    #########################################

    return list(ExpandedTorsion)


def IncludeH6(mol, ExpandedTorsion):
    for idx in ExpandedTorsion:
        atom = mol.GetAtomWithIdx(idx)
        for at in atom.GetNeighbors():
            if at.GetAtomicNum() == 1:
                ExpandedTorsion.append(at.GetIdx())

    ExpandedTorsion = set(ExpandedTorsion)
    return list(ExpandedTorsion)


def CapOpenValenced7(mol, ExpandedTorsion, TorsionQuartet):

    new_mol = Chem.RWMol()
    atom_map = {}
    new_atom_map = {}
    NOS_edged_idxs = {}
    C_edged_idxs = {}

    # add atoms
    atom_rank = 0
    quartet_new = []

    for idx in ExpandedTorsion:
        atom = mol.GetAtomWithIdx(idx)
        atom_map[idx] = new_mol.AddAtom(atom)
        if idx in TorsionQuartet:
            quartet_new.append(atom_rank) # transform idx from mol to new_mol with new idx
        atom_rank+=1

        count = 0
        neigh_idxs = findneighbour(mol, idx)
        for neigh_idx in neigh_idxs:
            if neigh_idx not in ExpandedTorsion and mol.GetAtomWithIdx(neigh_idx).GetAtomicNum() > 1:
                count += int(mol.GetBondBetweenAtoms(idx, neigh_idx).GetBondTypeAsDouble()) # single 1.0, double 2.0, aromatic 1.5, but with kekulize, aromatic will be either 1 or 2
        if (atom.GetAtomicNum() in [7, 8, 16]): #N,O,S
            NOS_edged_idxs[idx] = count
        if (atom.GetAtomicNum() in [6]): #C, C=C problem not considered
            C_edged_idxs[idx] = count

    # add bonds
    for idx in ExpandedTorsion:
        a = mol.GetAtomWithIdx(idx)
        for b in a.GetNeighbors():
            if b.GetIdx() not in ExpandedTorsion:
                continue
            bond = mol.GetBondBetweenAtoms(a.GetIdx(), b.GetIdx())
            bt = bond.GetBondType()
            if a.GetIdx() < b.GetIdx():  # Cauz each bond is enumerated twice
                new_mol.AddBond(atom_map[a.GetIdx()], atom_map[b.GetIdx()], bt)
    
    # cap H on C
    for idx,count in C_edged_idxs.items():
        for i in range(count):
            new_atom_map[idx * (i+1)*10] = new_mol.AddAtom(Chem.Atom(1))
            new_mol.AddBond(
                atom_map[idx], new_atom_map[idx * (i+1) * 10], Chem.BondType.SINGLE
            )

    # cap CH3 on N/O/S
    for idx,count in NOS_edged_idxs.items():
        for i in range(count):
            new_atom_map[idx * (i+1)*10] = new_mol.AddAtom(Chem.Atom(6))
            new_mol.AddBond(
                atom_map[idx], new_atom_map[idx * (i+1) * 10], Chem.BondType.SINGLE
            )
            for j in range(3): # add hydrogens in CH3
                new_atom_map[idx * (i+1)*10 + j +1] = new_mol.AddAtom(Chem.Atom(1))
                new_mol.AddBond(
                    new_atom_map[idx * (i+1)*10], new_atom_map[idx * (i+1) * 10 + j + 1], Chem.BondType.SINGLE
                )

    ## reorder quartet_new
    new_mol_mol = new_mol.GetMol()
    idxr1, idxr2, idxr3, idxr4 = checkbonded(new_mol_mol, quartet_new)
    quartet_new = [idxr1, idxr2, idxr3, idxr4]
    return new_mol.GetMol(), quartet_new

def GenStartingConf8(mol, quartet_new, outpath=Path("outputs"), name="test.sdf"):
    # Get initial 3D structure
    Chem.SanitizeMol(mol)
    Chem.Kekulize(mol)
    mol = Chem.AddHs(mol)
    
    quartet_new_str = ' '.join([str(idx) for idx in quartet_new])
    mol.SetProp("_Name", name)
    mol.SetProp("TORSION_ATOMS_FRAGMENT", quartet_new_str)

    # EmbedMolecule(mol) # initial single 3D structure
    # MMFFOptimizeMolecule(mol) # MMFF optimize

    params = ETKDGv3()
    params.numThreads=4
    params.useSmallRingTorsions=True
    params.pruneRmsThresh=0.1
    params.clearConfs=True
    
    cids = EmbedMultipleConfs(mol, numConfs=5, params=params) # 5 confs
    MMFFOptimizeMoleculeConfs(mol) # optimize, otherwise failed in xtb opt
    mol_confs = [ mol for cid in cids ]

    if not outpath.exists():
        outpath.mkdir()
    # writer = Chem.SDWriter(str(outpath / (name + ".sdf")))
    # writer.write(mol)
    # writer.close()

    writer = Chem.SDWriter(str(outpath / (name + ".sdf")))
    for i, mol_conf in enumerate(mol_confs):
        writer.write(mol_conf, confId=i)
    writer.close()

    return mol

class TorsionFragmentGenerator(object):
    def __init__(self, mol, outpath, name):
        rings = GetRingSystems(mol, includeSpiro=False)  # no spiro ring
        fcgps = [ mol.GetSubstructMatches(Chem.MolFromSmarts(fcgp)) for fcgp in RDK_NCOPS_group_SMARTS_NOS_simplified ]
        fcgp_flatten = []
        for fcgp in fcgps:
            for sub_fcgp in fcgp:
                sub_fcgp = sorted(sub_fcgp)
                if len(sub_fcgp)>1 and sub_fcgp not in fcgp_flatten:
                    fcgp_flatten.append(sub_fcgp)

        # mol_noH = copy.deepcopy(mol)
        # mol_noH = Chem.RemoveHs(mol_noH)  # RemoveHs at first
        mol_noH = Chem.RemoveHs(mol)

        TorsionPairs = GetTorsion0(mol_noH)
        print(">>> TorsionPairs:")
        print(TorsionPairs)
        out_names, new_mols = [], []

        print(">>> TorsionQuartet:")
        index = 0
        quartet_news, quartet_olds = [], []
        for TorsionPair in TorsionPairs:
            TorsionQuartet = GetQuartetAtoms1(mol_noH, TorsionPair)            
            if TorsionQuartet:
                print(TorsionQuartet)
                quartet_olds.append(TorsionQuartet) # Note that the molecule is canonicalized and thus the indexes may not be the same.
                ExpandedTorsion = ExpandQuartet2(mol_noH, TorsionQuartet)
                ExpandedTorsion = KeepRing3(mol_noH, ExpandedTorsion, rings)
                ExpandedTorsion = KeepOrtho4(mol_noH, ExpandedTorsion, TorsionQuartet, rings)
                ExpandedTorsion = KeepFuncGroup5(mol_noH, ExpandedTorsion, fcgp_flatten)
                ExpandedTorsion = IncludeH6(mol, ExpandedTorsion)
                # print(ExpandedTorsion)
                new_mol, quartet_new = CapOpenValenced7(mol, ExpandedTorsion, TorsionQuartet)
                new_mol = GenStartingConf8(new_mol, quartet_new, outpath, name + "_" + str(index)) # with 3D structure

                new_mols.append(new_mol)
                quartet_news.append(quartet_new)
                out_names.append(name + "_" + str(index))
                index += 1

        self.mols = new_mols
        self.old_quartets = quartet_olds
        self.new_quartets = quartet_news
        self.outnames = out_names
