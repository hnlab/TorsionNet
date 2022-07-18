import sys
import copy
import argparse
import re
import pickle
import operator
from pathlib import Path, PosixPath
import numpy as np
import scipy
from scipy import interpolate
from scipy.optimize import minimize
from matplotlib import pyplot as plt

from rdkit import Chem
from rdkit.Geometry import Point3D
from ase.units import Hartree, kcal, mol
kcalmol = kcal / mol
sys.path.append("/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/utils")
from GOutils import SPout2energy,GOout2energy

"""
To write out each mol's torsions' info into single .sdf file
"""

def xyzfile2coords(xyzfile):
    coords_text = xyzfile.read_text().split("\n")[2:-1]
    coords = [ coord_text.split()[1:] for coord_text in coords_text]
    new_atom_ps = [[float(one_coord) for one_coord in coord ] for coord in coords]

    return new_atom_ps

def SetNewPositions(mol,new_coords):
    """
    ::param1 mol: rdMol
    ::param2 coords: list of coords
    """
    # print(mol.GetConformer().GetAtomPosition(0).x) #raw
    mol_ = copy.deepcopy(mol) #not change raw mol
    conf = mol_.GetConformer()
    for i in range(mol_.GetNumAtoms()):
        x,y,z = new_coords[i]
        conf.SetAtomPosition(i,Point3D(x,y,z)) # will write into mols

    # print(mol.GetConformer().GetAtomPosition(0).x) #has changed
    return mol_

def curve_fit(angles,energies,subdir,outpath):
    func = interpolate.interp1d(x=angles, y=energies, kind="cubic")
    probe_angs = list(range(-180,181,1))
    pred_energies = [func(ang) for ang in probe_angs]
    min_E = min(pred_energies)
    min_ang = probe_angs[pred_energies.index(min_E)]
    print(min_ang)
    try:
        res = minimize(func,min_ang) # min_ang is the initial guess (as local minimum)
        # other options:
        # method: If not given, chosen to be one of BFGS, L-BFGS-B, SLSQP, depending on whether or not the problem has constraints or bounds.
    
        # But actually, because of the step of constrain but not fix in xtb optimization step,
        # There should be error which cannot be cancelled by minimization
        if res.success:
            min_ang = res.x[0] # min_ang
            min_E = res.fun if not isinstance(res.fun, list) else res.fun[0]# min_energy
    except ValueError: # in case out of range(-180, 180)
        min_ang = min_ang
        min_E = min_E

    rel_E = energies - min_E

    plt.figure(figsize=(14, 12), dpi=80)
    plt.plot(angles, rel_E, color="black", marker="o", label="calculation")
    plt.plot(probe_angs, np.array([ func(ang) for ang in probe_angs ]) - min_E, color = "blue", marker = "^", label="fitted curve")
    font1 = {"family": "Helvetica", "weight": "normal", "size": 30}
    font2 = {"family": "Helvetica", "weight": "normal", "size": 23}

    plt.title(subdir, font1)
    plt.legend(loc="upper left", prop=font2)
    plt.tick_params(labelsize=23)
    plt.xlabel("dihedral_angle", font2)
    plt.ylabel("$dE(kcal/mol)$", font2)
    plt.xlim((-180, 180))
    plt.savefig(outpath / f"{subdir}.png")
    plt.close()

    with open(outpath / f"{subdir}.func.pkl", 'wb') as f:
        pickle.dump(func, f)
    with open(outpath / f"{subdir}.min.csv", 'w') as f:
        f.write(f"{min_ang},{min_E}\n")


# if __name__ == "__main__":
if True:
    parser = argparse.ArgumentParser("write new coords to sdf file")
    parser.add_argument("--rootpath", type=PosixPath, help="absolute path for rootpath")
    parser.add_argument("--subdir", type=str, help="name of subdir")
    args = parser.parse_args()

    # rootpath = Path("/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/outputs/ShuoGu_D4_256")
    rootpath = args.rootpath
    startpath = rootpath / "0start"
    MMscanpath = rootpath / "1MMscans"
    xtbGOpath = rootpath / "2xtbGO"
    QMSPpath = rootpath / "3QMSP"
    QMGOpath = rootpath / "4QMGO"

    outmmsdfpath = rootpath / "outmmsdf"
    outqmsdfpath = rootpath / "outqmsdf"
    outqmfuncpath = rootpath / "outfunc"
    for path in [outmmsdfpath, outqmsdfpath, outqmfuncpath]:
        if not path.exists():
            path.mkdir()

    subdir = args.subdir
    MMscanpath_sub = MMscanpath / subdir
    xtbGOpath_sub = xtbGOpath / subdir
    QMGOpath_sub = QMGOpath / subdir
    QMSPpath_sub = QMSPpath / subdir
    outqmfuncpath_sub = outqmfuncpath / subdir
    if not outqmfuncpath_sub.exists():
        outqmfuncpath_sub.mkdir()

    print(f">>> Dealing with {subdir}")
    xtboptmols, QMoptmols = [], []
    onegroupsdf = MMscanpath_sub.glob("*-1_*.sdf")
    lines = (QMGOpath/"QMGO.csv").read_text().split("\n")
    wholegroupsdf = []
    for line in lines:
        if re.search(subdir, line):
            # name = line.split(",")[0].split(".opt.sdf.orcainp.xyz")[0]
            name = line.split(",")[0].split(".opt.sdf.orcainp.orcainp.xyz")[0]
            wholegroupsdf.append(MMscanpath_sub / name)

    """
    for sdffile in onegroupsdf:
    # for out
        mol = Chem.SDMolSupplier(str(sdffile), removeHs=False)[0]

        xtboptoutxyz = xtbGOpath_sub / (sdffile.name + ".opt.xyz")
        # QMoptoutxyz = QMGOpath_sub / (sdffile.name + ".opt.xyz.orcainp.xyz")

        ang = int(sdffile.name.split(".sdf")[0].split("_")[-1])
    
        xtbopt_coords = xyzfile2coords(xtboptoutxyz)
        # QMopt_coords = xyzfile2coords(QMoptoutxyz)

        xtbopt_mol = SetNewPositions(mol,xtbopt_coords)
        # QMopt_mol = SetNewPositions(mol,QMopt_coords)

        xtboptoutElog = QMSPpath_sub / (sdffile.name + ".opt.xyz.orcainp.log")
        # QMoptoutElog = QMGOpath_sub / (sdffile.name + ".opt.xyz.orcainp.log")

        xtbopt_E = SPout2energy(xtboptoutElog) * Hartree / kcalmol
        # QMopt_E = GOout2energy(QMoptoutElog) * Hartree / kcalmol

        xtbopt_mol.SetProp("_Name", sdffile.name + ".xtbopt")
        # QMopt_mol.SetProp("_Name", sdffile.name + ".QMopt")

        xtbopt_mol.SetProp("TORSION_ANGLE", str(ang))
        # QMopt_mol.SetProp("TORSION_ANGLE", str(ang))
        xtbopt_mol.SetProp("r2SCAN-3c_energy(kcal/mol)", "%.2f" % xtbopt_E)
        # QMopt_mol.SetProp("r2SCAN-3c_energy(kcal/mol)", "%.2f" % QMopt_E)

        xtboptmols.append(xtbopt_mol)
        # QMoptmols.append(QMopt_mol)

    writer = Chem.SDWriter(str(outmmsdfpath / f"{subdir}_one.sdf"))
    for ang in range(-180, 180, 15):
        for mol in xtboptmols:
            if mol.GetProp("TORSION_ANGLE") == str(ang):
                writer.write(mol)
    writer.close()

    # writer = Chem.SDWriter("test_QMopt_one.sdf")
    # for ang in range(-180, 180, 15):
        # for mol in QMoptmols:
            # if mol.GetProp("TORSION_ANGLE") == str(ang):
                # writer.write(mol)
    """
    
    angs, QMopt_Es = [], []
    if len(wholegroupsdf) == 24:
        for sdffile in wholegroupsdf:
            mol = Chem.SDMolSupplier(str(sdffile), removeHs=False)[0]

            # xtboptoutxyz = xtbGOpath_sub / (sdffile.name + ".opt.xyz")
            QMoptoutxyz = QMGOpath_sub / (sdffile.name + ".opt.xyz.orcainp.xyz")
            if not QMoptoutxyz.exists():
                QMoptoutxyz = QMGOpath_sub / (sdffile.name + ".opt.sdf.orcainp.xyz")

            ang = int(sdffile.name.split(".sdf")[0].split("_")[-1])
            angs.append(ang)
            # xtbopt_coords = xyzfile2coords(xtboptoutxyz)
            QMopt_coords = xyzfile2coords(QMoptoutxyz)

            # xtbopt_mol = SetNewPositions(mol,xtbopt_coords)
            QMopt_mol = SetNewPositions(mol,QMopt_coords)

            # xtboptoutElog = QMSPpath_sub / (sdffile.name + ".opt.xyz.orcainp.log")
            QMoptoutElog = QMGOpath_sub / (sdffile.name + ".opt.xyz.orcainp.log")
            if not QMoptoutElog.exists():
                QMoptoutElog = QMGOpath_sub / (sdffile.name + ".opt.sdf.orcainp.log")

            # xtbopt_E = SPout2energy(xtboptoutElog) * Hartree / kcalmol
            QMopt_E = GOout2energy(QMoptoutElog) * Hartree / kcalmol
            QMopt_Es.append(QMopt_E)

            # xtbopt_mol.SetProp("_Name", sdffile.name + ".xtbopt")
            QMopt_mol.SetProp("_Name", sdffile.name + ".QMopt")

            # xtbopt_mol.SetProp("TORSION_ANGLE", str(ang))
            QMopt_mol.SetProp("TORSION_ANGLE", str(ang))
            # xtbopt_mol.SetProp("r2SCAN-3c_energy(kcal/mol)", "%.2f" % xtbopt_E)
            QMopt_mol.SetProp("r2SCAN-3c_energy(kcal/mol)", "%.2f" % QMopt_E)

            # xtboptmols.append(xtbopt_mol)
            QMoptmols.append(QMopt_mol)

        # writer = Chem.SDWriter("test_xtbopt_whole.sdf")
        # for ang in range(-180, 180, 15):
            # for mol in xtboptmols:
                # if mol.GetProp("TORSION_ANGLE") == str(ang):
                    # writer.write(mol)
        # writer.close()

        writer = Chem.SDWriter(str(outqmsdfpath / f"{subdir}_whole.sdf"))
        for ang in range(-180, 180, 15):
            for mol in QMoptmols:
                if mol.GetProp("TORSION_ANGLE") == str(ang):
                    writer.write(mol)

        ang_E_dict = dict(zip(angs, QMopt_Es))
        sorted_by_ang = sorted(iter(ang_E_dict.items()), key=operator.itemgetter(0))
        ANG = np.array([q[0] for q in sorted_by_ang], dtype=int)
        ANG = np.append(ANG, 180)
        E = np.array([q[1] for q in sorted_by_ang], dtype=float)
        E = np.append(E, E[0])
        curve_fit(ANG, E, subdir, outqmfuncpath_sub)
    
    else:
        print(f"something wrong with {subdir}, take care!")

    # print(QMopt_Es)
    print(f">>> Finished with {subdir}")
    print()
    