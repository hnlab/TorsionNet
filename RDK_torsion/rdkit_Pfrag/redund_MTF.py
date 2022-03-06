import os
import copy
import pandas as pd
from pathlib import Path
from collections import defaultdict
from rdkit import Chem


TORSIONCHECKER="/pubhome/soft/ZBH/Torsion Analyzer/torsionchecker"
tor_lic = Path("/pubhome/qcxia02/0scripts/torsionchecker_lic.txt").read_text().split("\n")[0]
tor_lib = "/pubhome/soft/ZBH/Torsion Analyzer/tor_lib.xml" # tor_lib_v2

"""
This script is now specific to AmpC subset compounds
"""

def tor_lib_SMARTS(sdffile,smi,outpath):
    # note that is recommended that sdffile contains only on molecule
    mol = Chem.SDMolSupplier(str(sdffile), removeHs = False)[0]
    raw_mol = copy.deepcopy(mol)

    # Not necessary
    # canonicalize
    mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol))
    writer = Chem.SDWriter("tmp.sdf")
    writer.write(mol)
    writer.close()

    # Use torsionchecker to obtain torsion SMARTS for pattern mapping #
    # Can be substituted by any other tools that can obtain torsion SMARTS"
    # os.system(f'\"{TORSIONCHECKER}\" -m {str(sdffile)} -t \"{tor_lib}\" -r tmp.csv -i {tor_lic}')
    os.system(f'\"{TORSIONCHECKER}\" -m tmp.sdf -t \"{tor_lib}\" -r tmp.csv -i {tor_lic}')

    df = pd.read_csv("tmp.csv", sep="\t")
    tor_smarts = list(df["SMARTS"])
    torsion_matches = [ raw_mol.GetSubstructMatches(Chem.MolFromSmarts(tor_smart))[0] for tor_smart in tor_smarts ] # for multiple matches, only take the first one
    torsion_frag = set([ str(idx) for idx in raw_mol.GetProp("TORSION_ATOMS_FRAGMENT").split() ])
    torsion_matches = [ set(match) for match in torsion_matches]

    mapped_smarts = ""
    for i,match in enumerate(torsion_matches):
        if torsion_frag.issubset(match):
            mapped_smarts = tor_smarts[i] 
            break # if map one, then finish

    os.system("rm tmp.sdf tmp.csv")
    # os.system("rm tmp.csv")

    ###################################################################
    if mapped_smarts:
        print(mapped_smarts)
        raw_mol.SetProp("SMILES",smi) # save things
        raw_mol.SetProp("Torsion_SMARTS", mapped_smarts)
    
        writer = Chem.SDWriter(str(outpath / sdffile.name))
        writer.writer(mol)
        writer.close()
    return

if True:
    sdfpath = Path("/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/outputs")
    
    smilist = []
    filt_files = []
    
    for i in range(59):
        files = list(sdfpath.glob(f"cpds.binders.{i}_*.sdf"))
        filenames = [ file.name for file in files]
        smis = [ Chem.MolToSmiles(Chem.SDMolSupplier(str(file))[0]) for file in files ]
        smi_filename_dict = defaultdict(list)
        # for smi, filename in zip(smis, filenames):
            # smi_filename_dict[smi].append(filename)
        for smi, sdffile in zip(smis, files):
            smi_filename_dict[smi].append(sdffile)
        for smi in smis:
            if smi not in smilist:
                smilist.append(smi)
                # filt_files.append(smi_filename_dict[smi])
                filt_files.extend(smi_filename_dict[smi])

    
    for i in range(200):
        files  = list(sdfpath.glob(f"cpds.nonbinders200.{i}_*.sdf"))
        filenames = [ file.name for file in files]
        smis = [ Chem.MolToSmiles(Chem.SDMolSupplier(str(file))[0]) for file in files ]
        smi_filename_dict = defaultdict(list)
        # for smi, filename in zip(smis, filenames):
            # smi_filename_dict[smi].append(filename)
        for smi, sdffile in zip(smis, files):
            smi_filename_dict[smi].append(sdffile) # ensure different smarts in same smi are saved

        for smi in smis:
            if smi not in smilist:
                for sdffile in smi_filename_dict[smi]:
                    tor_lib_SMARTS(sdffile, smi, outpath="redund")

                smilist.append(smi)
                # filt_files.append(smi_filename_dict[smi])
                filt_files.extend(smi_filename_dict[smi])

    
    print(len(filt_files))
