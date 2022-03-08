#!/pubhome/qcxia02/miniconda3/envs/basic/bin/python3.9
import os
import sys
import argparse
from pathlib import Path, PosixPath

from rdkit import Chem
from rdkit.Chem import Draw

sys.path.append("/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/utils")
from TFG import TorsionFragmentGenerator


if True:
    """
    Usage: python rdkit_Pfrag.py --mol2 mol2file.mol2 --outpath outpath --imgpath imgpath
    """
    parser = argparse.ArgumentParser("RDKit Torsion Fragment Generator")
    parser.add_argument("--sdf", type=str, help="path for sdf file with single M")
    parser.add_argument("--mol2", type=str, help="path for mol2 file with single M")
    parser.add_argument("--outpath",type=PosixPath, help="path for output .sdf mols")
    parser.add_argument("--imgpath",type=PosixPath, help="path for output imgs")
    args = parser.parse_args()

    # mol = Chem.MolFromSmiles("COc3ccc(C)c(C(C)NC(=O)Nc2nc(c1ccncc1)cs2)c3C")
    # mol = Chem.SDMolSupplier(args.sdf, removeHs = False)[0]
    # mol2file = "/pubhome/qcxia02/git-repo/DL4molcst/scripts/dock_models/3-conf-docks/results/AmpC/compounds/mol2mols/mol2mols/cpds.binders.0.101.mol2"
    mol2file = args.mol2
    sdffile = args.sdf
    outpath = args.outpath
    imgpath = args.imgpath

    for path in [outpath, imgpath]:
        if not path.exists():
            os.system(f"mkdir -p {str(path)}")

    if args.sdf:
        mol = Chem.SDMolSupplier(mol2file, removeHs = False)[0] # readfile
        molname = Path(sdffile).name[:-4]

    if args.mol2:
        mol = Chem.MolFromMol2File(mol2file, removeHs = False) # readfile
        molname = Path(mol2file).name[:-5]

    print(f">>>>>> Dealing with {molname} <<<<<<")

    mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol)) # canonicalize
    Chem.SanitizeMol(mol)
    Chem.Kekulize(mol, clearAromaticFlags=True) # To bypass complex Kekulize Error
    mol = Chem.AddHs(mol)

    ### OUTPUT .sdf molecules will be saved in outpath dir at current working directory
    TF = TorsionFragmentGenerator(mol, outpath=outpath, name=molname) # This name is the prefix for output SDF mol
    # TF = TorsionFragmentGenerator(mol, outpath="outputs", name="test") # This name is the prefix for output SDF mol


    new_mols = TF.mols
    old_quartets = TF.old_quartets
    new_quartets = TF.new_quartets
    out_names = TF.outnames

    print(f">>> Total {len(new_mols)} torsion fragments")

    for i, new_mol in enumerate(new_mols):
        quartet = new_mol.GetProp("TORSION_ATOMS_FRAGMENT").split(" ")
        quartet = [int(idx) for idx in quartet]
        print(f">>>>>> New quartet atoms for {i}:")
        print(quartet)
    
    if new_quartets: # only when with torsion, can output
        img = Draw.MolsToGridImage(
        [ Chem.MolFromSmiles(Chem.MolToSmiles(new_mol)) for new_mol in new_mols],
        molsPerRow=5,
        subImgSize=(350,350),
        legends = out_names
    )
        img.save(str(imgpath / (molname + ".png")))

    print(f">>>>>> Finished with {molname} <<<<<<")
    print()
