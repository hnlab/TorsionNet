import sys
import copy
import pandas as pd
import pickle
import argparse
from pathlib import Path
from rdkit import Chem
from rdkit.Chem.Lipinski import RotatableBondSmarts
from rdkit.Chem.rdMolTransforms import GetDihedralDeg,SetDihedralDeg

sys.path.append("/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/utils")
# from utils import GetRingSystems, findneighbour, Get_sorted_heavy
from TFG import GetTorsion0, GetQuartetAtoms1

def getdihedralangle(mol, torsionquartet):
    idxr1, idxr2, idxr3, idxr4 = torsionquartet
    torsionangle = GetDihedralDeg(mol.GetConformer(), idxr1, idxr2, idxr3, idxr4,)
    return torsionangle

def calc_torsion(mol, torsionquartet, func, minE):
    torsionangle = getdihedralangle(mol, torsionquartet)
    torsion_strain = func(torsionangle) - minE
    return torsion_strain

def write_best_mol(mol, torsionquartets, torsionangles, outmolpath):
    for i in range(len(torsionquartets)):
        torsionquartet = torsionquartets[i]
        idxr1, idxr2, idxr3, idxr4 = torsionquartet
        torsionangle = torsionangles[i]
        SetDihedralDeg(mol.GetConformer(), idxr1, idxr2, idxr3, idxr4, torsionangle)
        mol.SetProp(f"Torsion_Quartet{i+1}",f"{idxr1} {idxr2} {idxr3} {idxr4}")
        mol.SetProp(f"Torsion_angles{i+1}", "%.2f" % torsionangle)
    writer = Chem.SDWriter(str(outmolpath / "test.sdf"))
    writer.write(mol)
    writer.close()
    return

def removeNone(probe_list):
    while None in probe_list:
        probe_list.remove(None)
    return probe_list

if __name__ == "__main__":
    # parser = argparse.ArgumentParser()
    outfuncpath = Path("/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/outputs/ShuoGu_AmpC_44/outfunc")
    outmolpath = Path("/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/outputs/ShuoGu_AmpC_44/outbestmol")

    inp = "/pubhome/qcxia02/Downloads/dataset/ShuoGu/AmpC/ZINC000456896024.mol2"
    prefix = Path(inp).stem
    targetpath = list(outfuncpath.glob(f"{prefix}*"))
    # print(targetpath)
    funcs, minangles, minEs = [],[],[]
    for i in range(len(targetpath)):
        subname = prefix + f"_{i}"
        with open(outfuncpath / subname / f"{subname}.func.pkl", 'rb') as f:
            func = pickle.load(f)
        with open(outfuncpath / subname / f"{subname}.min.csv", 'r') as f:
            lines = f.readlines()
        minangle = float(lines[0].split(",")[0])
        minE = float(lines[0].split(",")[1])

        funcs.append(func)
        minangles.append(minangle)
        minEs.append(minE)
    mol = Chem.MolFromMol2File(inp, removeHs=False)
    canonical_mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol))
    Chem.SanitizeMol(canonical_mol)
    Chem.Kekulize(canonical_mol, clearAromaticFlags=True)

    raw_order = list(map(int, mol.GetProp("_smilesAtomOutputOrder")[1:-2].split(",")))
    for i,atom in enumerate(mol.GetAtoms()):
        if atom.GetAtomicNum() == 1: # remove H in raw_order through this way to keep mapping
            raw_order.remove(i)

    new_order = list(range(len(canonical_mol.GetAtoms())))
    if len(raw_order) == len(new_order):
        new_raw_mapping = dict(zip(new_order, raw_order))

        torsions = GetTorsion0(canonical_mol)
        TorsionQuartets_ = [ GetQuartetAtoms1(canonical_mol, torsion) for torsion in torsions ]
        TorsionQuartets =  removeNone(TorsionQuartets_)
        print(TorsionQuartets)
        raw_TorsionQuartets = [[ new_raw_mapping[index] for index in TorsionQuartet ] for TorsionQuartet in TorsionQuartets]
        print(raw_TorsionQuartets)
        TorsionAngles = [getdihedralangle(mol, raw_TorsionQuartets[i]) for i in range(len(raw_TorsionQuartets))]
        print(TorsionAngles)
        TorsionStrains = [calc_torsion(mol, raw_TorsionQuartets[i], funcs[i], minEs[i]) for i in range(len(raw_TorsionQuartets))]
        print(sum(TorsionStrains))
        print(max(TorsionStrains))

        write_best_mol(mol, raw_TorsionQuartets, minangles, outmolpath)
    else:
        print("something wrong! please check!")