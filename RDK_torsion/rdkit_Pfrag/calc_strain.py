"""
Note that TFG_TL.py does not need canonicalization, so I turned them off - 2022/05/27
"""
import sys
import copy
import pandas as pd
import pickle
import argparse
from pathlib import Path, PosixPath
from rdkit import Chem
from rdkit.Chem.Lipinski import RotatableBondSmarts
from rdkit.Chem.rdMolTransforms import GetDihedralDeg,SetDihedralDeg

sys.path.append("/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/utils")
# from utils import GetRingSystems, findneighbour, Get_sorted_heavy
# from TFG import GetTorsion0, GetQuartetAtoms1
# from TFG_TL import GetTorsionQuartet01
from TFG_TEU import GetTorsionQuartet01

def getdihedralangle(mol, torsionquartet):
    idxr1, idxr2, idxr3, idxr4 = torsionquartet
    torsionangle = GetDihedralDeg(mol.GetConformer(), idxr1, idxr2, idxr3, idxr4,)
    return torsionangle

def calc_torsion(mol, torsionquartet, func, minE):
    torsionangle = getdihedralangle(mol, torsionquartet)
    torsion_strain = func(torsionangle) - minE
    return torsion_strain

def write_best_mol(mol_, torsionquartets, torsionangles, outmolpath, molname):
    mol = copy.deepcopy(mol_)
    for i in range(len(torsionquartets)):
        torsionquartet = torsionquartets[i]
        idxr1, idxr2, idxr3, idxr4 = torsionquartet
        torsionangle = torsionangles[i]
        SetDihedralDeg(mol.GetConformer(), idxr1, idxr2, idxr3, idxr4, torsionangle)
        mol.SetProp(f"Torsion_Quartet{i+1}",f"{idxr1} {idxr2} {idxr3} {idxr4}")
        mol.SetProp(f"Min_Torsion_Angle{i+1}", "%.2f" % torsionangle)
    writer = Chem.SDWriter(str(outmolpath / (molname + "_best.sdf")))
    writer.write(mol)
    writer.close()
    return

def write_probe_mol(mol_, torsionquartets, torsionangles, outmolpath, molname, TorsionStrains):
    mol = copy.deepcopy(mol_)
    for i in range(len(torsionquartets)):
        torsionquartet = torsionquartets[i]
        idxr1, idxr2, idxr3, idxr4 = torsionquartet
        torsionangle = torsionangles[i]
        mol.SetProp(f"Torsion_Quartet{i+1}",f"{idxr1} {idxr2} {idxr3} {idxr4}")
        mol.SetProp(f"Torsion_Angle{i+1}", "%.2f" % torsionangle)
        mol.SetProp(f"Torsion_Strain{i+1}", "%.2f" % TorsionStrains[i])

    mol.SetProp("Sum_Torsion_Strain", "%.2f" % sum(TorsionStrains))
    mol.SetProp("Max_Torsion_Strain", "%.2f" % max(TorsionStrains))
    writer = Chem.SDWriter(str(outmolpath / (molname + "_probe.sdf")))
    writer.write(mol)
    writer.close()
    return

def removeNone(probe_list):
    while None in probe_list:
        probe_list.remove(None)
    return probe_list
"""
Usage: python calc_strain.py --rootpath $rootpath --mol2 $mol2file
"""
if __name__ == "__main__":
    parser = argparse.ArgumentParser("calc_strain")
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('--mol2', type=str, help="Absolute path of Mol2 file to gererate conformers for")
    input_group.add_argument('--sdf', type=str, help="Absolute path of SDF file to gererate conformers for")
    input_group.add_argument('--smiles', type=str, help="SMILES string of molecule")
    parser.add_argument("--rootpath", type=PosixPath, help="absolute path of rootpath of study", required=True)
    args = parser.parse_args()

    rootpath = args.rootpath
    outfuncpath = rootpath / "outfunc"
    outbestmolpath = rootpath / "outbestmol"
    outprobemolpath = rootpath / "outprobemol"
    # outfuncpath = Path("/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/outputs/ShuoGu_AmpC_44/outfunc")
    # outbestmolpath = Path("/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/outputs/ShuoGu_AmpC_44/outbestmol")
    # outprobemolpath = Path("/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/outputs/ShuoGu_AmpC_44/outprobemol")
    for path in [outbestmolpath, outprobemolpath]:
        if not path.exists():
            path.mkdir()
    
    ## Deal with mol like in TFG.py ##
    if hasattr(args, 'mol2') and args.mol2 is not None:
        inp = args.mol2
        prefix = Path(inp).stem
        mol = Chem.MolFromMol2File(inp, removeHs=False) # only read the first mol2 mol
        canonical_mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol))

    if hasattr(args, 'sdf') and args.sdf is not None:
        inp = args.sdf
        prefix = Path(inp).stem
        mol = Chem.SDMolSupplier(inp, removeHs=False)[0] # only read the first sdf mol
        canonical_mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol))

    # if hasattr(args, 'smiles') and args.mol2 is not None:
    #     prefix = args.smiles

    print(f">>> Dealing with {prefix}")
    
    Chem.SanitizeMol(canonical_mol)
    Chem.Kekulize(canonical_mol, clearAromaticFlags=True)

    raw_order = list(map(int, mol.GetProp("_smilesAtomOutputOrder")[1:-2].split(",")))
    for i,atom in enumerate(mol.GetAtoms()):
        if atom.GetAtomicNum() == 1: # remove H in raw_order through this way to keep mapping
            raw_order.remove(i)
    ##########################

    ## read func and energy ##
    # targetpath = list(outfuncpath.glob(f"{prefix}*"))
    targetpath = list(outfuncpath.glob(f"{prefix}_*")) #  to avoid extra mismatch
    # print(targetpath)
    funcdict = {}
    for subname in targetpath:
        index = int(subname.name.split("_")[-1])
        subname = subname.name
        with open(outfuncpath / subname / f"{subname}.func.pkl", 'rb') as f:
            func = pickle.load(f)
        with open(outfuncpath / subname / f"{subname}.min.csv", 'r') as f:
            lines = f.readlines()
        minangle = float(lines[0].strip().split(",")[0])
        minE = float(lines[0].strip().split(",")[1])

        funcdict[index] = {
            "func": func,
            "minangle": minangle,
            "minE": minE
        }
        TorsionStrains = []
    ##########################

    new_order = list(range(len(canonical_mol.GetAtoms())))
    if len(raw_order) == len(new_order):
    # if True:
        new_raw_mapping = dict(zip(new_order, raw_order))

        # torsions = GetTorsion0(canonical_mol)
        # TorsionQuartets_ = [ GetQuartetAtoms1(canonical_mol, torsion) for torsion in torsions ]
        # _, _, TorsionQuartets_ = GetTorsionQuartet01(mol)
        TorsionQuartets_, _ = GetTorsionQuartet01(mol)
        TorsionQuartets =  removeNone(TorsionQuartets_)
        print(TorsionQuartets)
        raw_TorsionQuartets = TorsionQuartets
        # raw_TorsionQuartets = [[ new_raw_mapping[index] for index in TorsionQuartet ] for TorsionQuartet in TorsionQuartets]
        # print(raw_TorsionQuartets)
        TorsionAngles = [getdihedralangle(mol, raw_TorsionQuartets[i]) for i in range(len(raw_TorsionQuartets))]
        print(TorsionAngles)
        # TorsionStrains = [calc_torsion(mol, raw_TorsionQuartets[i], funcs[i], minEs[i]) for i in range(len(raw_TorsionQuartets))]
        for i in range(len(raw_TorsionQuartets)):
            try:
                TorsionStrain = calc_torsion(mol, raw_TorsionQuartets[i], funcdict[i]["func"], funcdict[i]["minE"])
                TorsionStrains.append(TorsionStrain)
            except KeyError:
                TorsionStrains.append(0)
        print(TorsionStrains)

        with open(args.rootpath / "summary.csv",'a') as f:
            f.write(f"{prefix}\t{'%.2f' % sum(TorsionStrains)}\t{'%.2f' % max(TorsionStrains)}\n")
        print(sum(TorsionStrains))
        print(max(TorsionStrains))
        mol_ = copy.deepcopy(mol)
        write_probe_mol(mol_, raw_TorsionQuartets, TorsionAngles, outprobemolpath, prefix, TorsionStrains)
        # write_best_mol(mol_, raw_TorsionQuartets, minangles, outbestmolpath, prefix)

    else:
        print("something wrong! please check!")
    
    print(f">>> Finished with {prefix}")
    print("")
