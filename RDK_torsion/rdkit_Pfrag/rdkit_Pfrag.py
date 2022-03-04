#!/pubhome/qcxia02/miniconda3/envs/basic/bin/python3.9
import sys
import argparse
from rdkit import Chem

sys.path.append("/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/utils")
from TFG import TorsionFragmentGenerator


if True:
    parser = argparse.ArgumentParser("RDKit Torsion Fragment Generator")
    parser.add_argument("--sdf", type=str, help="absolute path for sdf file with single M")
    args = parser.parse_args()

    # mol = Chem.MolFromSmiles("COc3ccc(C)c(C(C)NC(=O)Nc2nc(c1ccncc1)cs2)c3C")
    mol = Chem.SDMolSupplier(args.sdf, removeHs = False)[0]
    Chem.SanitizeMol(mol)
    mol = Chem.AddHs(mol)

    TF = TorsionFragmentGenerator(mol, name="test") # This name is the prefix for output SDF mol

    new_mols = TF.mols
    new_quartets = TF.quartets
    out_names = TF.outnames

    print(f">>> Total {len(new_mols)} torsion fragments")

    for i, new_mol in enumerate(new_mols):
        quartet = new_mol.GetProp("TORSION_ATOMS_FRAGMENT").split(" ")
        quartet = [int(idx) for idx in quartet]
        print(f">>>>>> New quartet atoms for {i}:")
        print(quartet)
