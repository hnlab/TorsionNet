import os
import json
from pathlib import Path
from functools import partial
import xml.etree.ElementTree as ET #For reading and writing XML files
from matplotlib import pyplot as plt
from matplotlib.pyplot import MultipleLocator
import numpy as np
import pandas as pd
from scipy import interpolate
from scipy.optimize import minimize
from rdkit import Chem
from rdkit.Chem import rdMolTransforms
from ase.units import Hartree, kcal, mol
kcalmol = kcal / mol

def parseTorsionfromcsvline(csvlinelist):
    TS, TQ = [], []
    if len(csvlinelist) % 7 != 0:
        return np.nan, np.nan
    else:
        torsioncount = int(len(csvlinelist) / 7) - 1
        for i in range(torsioncount):
            if np.isnan(csvlinelist[i*7+7]) or i*7+7 == len(csvlinelist):
                # return TS, TQ
                break
            else:
                TS.append(csvlinelist[i*7+7+6])
                TQ.append(json.loads(csvlinelist[i*7+7+4]))
    return TS, TQ

def maketorsiondict(csvfile):
    df = pd.read_csv(csvfile, header=None)
    pdbid_torsion_dict = {}
    for index in range(len(df)): # ith row
        pdbid = df[0][index]
        TS, TQ = parseTorsionfromcsvline(list(df.iloc[index]))
        pdbid_torsion_dict[pdbid] = [TS, TQ]
    return pdbid_torsion_dict

def parseQMoptE(orcaoutfile):
    lines = Path(orcaoutfile).read_text().split("\n")
    energy = float(np.nan)
    for i,line in enumerate(lines):
        if "*** OPTIMIZATION RUN DONE ***" in line:
            energy = float(lines[i-3].split()[4])
    QMopt_E = energy * Hartree / kcalmol

    return QMopt_E

def curve_fit(angles,energies):
    func = interpolate.interp1d(x=angles, y=energies, kind="cubic")
    probe_angs = list(range(-180,181,1))
    pred_energies = [func(ang) for ang in probe_angs]
    min_E = min(pred_energies)
    min_ang = probe_angs[pred_energies.index(min_E)]
    # print(min_ang)
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

    return rel_E, func, min_E

def QMoptE_list(optlogpath,angles=list(range(-180,181,15))):
    print(optlogpath)
    QMoptEs = []
    # try:
    for ang in range(-180,180,15):
        strang = "%03d" % ang
        logfile = list(optlogpath.glob(f"*_{strang}.*.log"))[0]
        QMoptEs.append(parseQMoptE(logfile))
    QMoptEs.append(QMoptEs[0])
    # except:
        # print("qmopt parse error")
    rel_E, func, min_E = curve_fit(angles,np.array(QMoptEs))

    return rel_E, func, min_E

def CSD_hist(TorsionSmartsPattern, TorLibXMLfile):
    tree = ET.parse(TorLibXMLfile)
    root = tree.getroot()
    for torsionRule in root.iter("torsionRule"): #Loop again
        if TorsionSmartsPattern == torsionRule.get("smarts"):
            # bins = list(torsionRule.find("histogram_shifted").findall("bin")) # for analysis
            bins = list(torsionRule.find("histogram").findall("bin")) # for view
            bincounts = [ int(bin.get("count")) for bin in bins ]
            break
            
    return bincounts

def CSDTEU_list(TorsionSmartsPattern, TorLibXMLfile):
    print(TorsionSmartsPattern)
    from math import ceil
    tree = ET.parse(TorLibXMLfile)
    root = tree.getroot()
    for torsionRule in root.iter("torsionRule"): #Loop again
        if TorsionSmartsPattern == torsionRule.get("smarts"):
            # bins = list(torsionRule.find("histogram_shifted").findall("bin")) # for analysis
            bins = list(torsionRule.find("histogram_converted").findall("bin")) # for view
            hist_E = [ float(bin.get("energy")) for bin in bins ] # only for those "specific" and "exact"
            break
    probe_angs = list(range(-180,181,1))
    TEU_energies = []
    for theta in probe_angs:
        bin_num = ceil(theta / 10) + 17
        energy = (hist_E[bin_num]-hist_E[(bin_num+35)%36])/10.0*(theta-(bin_num-17)*10)+hist_E[bin_num]
        TEU_energies.append(energy)

    return TEU_energies

def sdf_ang_list(sdffile, torquartet):
    mols = Chem.SDMolSupplier(sdffile)
    idxr1,idxr2,idxr3,idxr4 = torquartet
    degs = []
    for mol in mols:
        deg = rdMolTransforms.GetDihedralDeg(mol.GetConformer(),idxr1,idxr2,idxr3,idxr4)
        degs.append(deg)
    return degs

def get_map_index(refsdffile, probesdffile):
    mol1 = Chem.SDMolSupplier(refsdffile)[0]
    mol2 = Chem.SDMolSupplier(probesdffile)[0]
    
    canonical_mol1 = Chem.MolFromSmiles(Chem.MolToSmiles(mol1,isomericSmiles=False))
    # print(Chem.MolToSmiles(canonical_mol1))
    # print(list(map(int, mol1.GetProp("_smilesAtomOutputOrder")[1:-2].split(","))))
    
    canonical_mol2 = Chem.MolFromSmiles(Chem.MolToSmiles(mol2,isomericSmiles=False))
    # print(Chem.MolToSmiles(canonical_mol2))
    # print(list(map(int, mol2.GetProp("_smilesAtomOutputOrder")[1:-2].split(","))))
    
    mapindex = dict(zip(
        list(map(int, mol1.GetProp("_smilesAtomOutputOrder")[1:-2].split(","))),
        list(map(int, mol2.GetProp("_smilesAtomOutputOrder")[1:-2].split(","))),
        ))
    return mapindex
