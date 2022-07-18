#!/pubhome/qcxia02/miniconda3/envs/basic/bin/python3.9
from rdkit import Chem
from rdkit.Chem.Lipinski import RotatableBondSmarts
from rdkit.Chem.rdDistGeom import EmbedMolecule, EmbedMultipleConfs, ETKDGv3
from rdkit.Chem.rdForceFieldHelpers import MMFFOptimizeMolecule,MMFFOptimizeMoleculeConfs
import sys
from pathlib import Path

def RDK_sample(mol, rmsd=0.1, numConfs=20, outpath=Path("."), name="test"):
    params = ETKDGv3()
    params.numThreads=4
    params.useSmallRingTorsions=True
    params.pruneRmsThresh=rmsd
    params.clearConfs=True
    
    cids = EmbedMultipleConfs(mol, numConfs=numConfs, params=params) # 5 confs
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

if __name__ == "__main__":
    sdffile = sys.argv[1]
    print(sdffile)
    name = sys.argv[2]
    mol = Chem.SDMolSupplier(str(sdffile), removeHs=False)[0]

    rmsd=0.1
    numConfs=20
    outpath = Path(".")

    RDK_sample(mol,rmsd,numConfs,outpath,name)
