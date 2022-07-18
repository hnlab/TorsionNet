from pathlib import Path, PosixPath
import argparse

from rdkit import Chem
from rdkit.Chem import rdMolTransforms
from rdkit.Chem import rdmolops

def findneighbour(mol, idxr) -> list:
    atom = mol.GetAtomWithIdx(idxr)  # Atoms: c(harged)atom
    neighbours = atom.GetNeighbors()  # list of Atoms: neighbour atoms of catom
    neigh_idxr = [
        x.GetIdx() for x in neighbours
    ]  # list of int: idxes of neighbours, include hydrogens
    # neigh_idxr = [x.GetIdx() for x in neighbours if x.GetAtomicNum != 1] #list of int: idxes of neighbours, without hydrogens
    return neigh_idxr

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

def clash(mol, torsion_quartet, cutoff):
    idxr1, idxr2, idxr3, idxr4 = checkbonded(mol, torsion_quartet)
    torsion_quartet_new = f"{idxr1} {idxr2} {idxr3} {idxr4}"
    mol.SetProp("TORSION_ATOMS_FRAGMENT", torsion_quartet_new)
    clash = False
    for angle in range(-180, 180, 15):
        rdMolTransforms.SetDihedralDeg(mol.GetConformer(),idxr1,idxr2,idxr3,idxr4,angle)
        dmat = rdmolops.Get3DDistanceMatrix(mol)
        for item in dmat[dmat<cutoff]: # cutoff set to <1A here
            if item!=0:
                clash = True
    return clash

def RDKit_rotate(mol, order, outpath, torsion_quartet):
    idxr1, idxr2, idxr3, idxr4 = checkbonded(mol, torsion_quartet)
    # torsion_quartet_new = f"{idxr1} {idxr2} {idxr3} {idxr4}"
    # mol.SetProp("TORSION_ATOMS_FRAGMENT", torsion_quartet_new)

    outpath = Path(outpath) / mol.GetProp("_Name")
    if not outpath.exists():
        outpath.mkdir()

    for ang in range(-180, 180, 15):
        rdMolTransforms.SetDihedralDeg(
            mol.GetConformer(), int(idxr1), int(idxr2), int(idxr3), int(idxr4), ang
        )

        w = Chem.SDWriter(
            str(outpath / f'{mol.GetProp("_Name")}-{order+1}_{"%03d" % ang}.sdf')
        )
        w.write(mol)


if True:
    parser = argparse.ArgumentParser()
    parser.add_argument("--inpath", type=PosixPath, help="absolute path of input fragments")
    parser.add_argument("--outpath", type=PosixPath, help="absolute path of output fragments")
    parser.add_argument("--torsionquartet", nargs="+", help="set torsion quartet manually")
    args = parser.parse_args()
    inpath = args.inpath
    outpath = args.outpath
    # inpath = Path("/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/outputs/0start")
    # outpath = Path("/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/outputs/1MMscans")
    if not outpath.exists():
        outpath.mkdir()
    
    for sdffile in inpath.iterdir():
        if sdffile.name.endswith(".sdf"):
            mols = Chem.SDMolSupplier(str(sdffile), removeHs=False)
            for i, mol in enumerate(mols):
                try:
                    torsion_quartet = mol.GetProp("TORSION_ATOMS_FRAGMENT").split()
                except:
                    if args.torsionquartet:
                        torsion_quartet = [int(i) for i in args.torsionquartet]
                print(torsion_quartet)
                # if not clash(mol, torsion_quartet, cutoff=0.9): # not used because I want more reasonable structures at all angles
                RDKit_rotate(mol, order = i, outpath = outpath, torsion_quartet = torsion_quartet)
