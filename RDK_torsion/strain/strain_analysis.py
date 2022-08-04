import os 
import sys
from multiprocessing import Pool
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

sys.path.append("/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/strain/utils")
from utils import QMoptE_list, CSD_hist, CSDTEU_list, sdf_ang_list, get_map_index, maketorsiondict
from plot import plotall

def plotallall(smarts,bincounts,func,min_E,rel_E,TEU_energies,degstatistics,xtaldeg,bestnofilter_deg,bestfilter_deg,topnofilter_deg,topfilter_deg,cate,imgpath,rowidx,axs):
    for num in range(5):
        if num == 0:
            # 1. plot_TEU
            ax1 = axs[rowidx][num]
            ax1.set_ylim(0,10000)
            x = list(range(-175,180,10))
            ax1.bar(x,bincounts, width=10, color="b",align="center",edgecolor="black",label="CSD histogram")
            # 2. plot_QMoptE
            angles = list(range(-180,181,15))
            probe_angs = list(range(-180,181,1))
            ax2 = ax1.twinx()
            lns1 = ax2.plot(probe_angs, np.array([ func(ang) for ang in probe_angs ]) - min_E, color = "black",label="QM fragment", linewidth=5)
            ax2.plot(angles, rel_E, color="red", marker="^", markersize=10,linestyle="")
            ax2.set_ylim(0,20)

            # 3. plot CSD-TEU
            lns2 = ax2.plot(probe_angs, TEU_energies, color="blue", label="CSD TEU",linewidth=5)
            lns = lns1+lns2
            labs = [l.get_label() for l in lns]

            # plt.close()
            plt.xlim(-185,185)
            x_major_locator = MultipleLocator(30) # The sep is set to 30
            ax = plt.gca() # the instance of 2 axes
            ax.xaxis.set_major_locator(x_major_locator)
            # plt.figure()
            font1 = {"family": "Helvetica", "weight": "normal", "size": 30}
            font2 = {"family": "Helvetica", "weight": "normal", "size": 23}
            ax1.set_title(smarts, font1)
            ax1.set_xlabel("Angle(degree)",font2)
            ax1.set_ylabel("CSD histogram count",font2)
            ax2.set_ylabel("QM relE(kcal/mol) or CSD TEU",font2)
            # ax2.plot(probe_angs, TEU_energies, color="black", marker="o", markersize=5,linestyle="")
            ax1.tick_params(labelsize=23)
            ax2.tick_params(labelsize=23)
            ax2.legend(lns, labs, loc="upper right",prop={"size":23})
            
        else:
            ax1 = axs[rowidx][num]
            degs = degstatistics[num-1]
            ax1.hist(degs, bins=int(360/10), range=[-180,180],edgecolor="black")
            ax1.set_ylim(0,100)
            ax2 = ax1.twinx()
            lns1 = ax2.plot(probe_angs, np.array([ func(ang) for ang in probe_angs ]) - min_E, color = "black",label="QM fragment", linewidth=5)
            ax2.plot(angles, rel_E, color="red", marker="^", markersize=10,linestyle="")
            ax2.set_ylim(0,20)
            lns2 = ax2.plot(probe_angs, TEU_energies, color="blue", label="CSD TEU",linewidth=5)
            ax1.tick_params(labelsize=23)
            ax2.tick_params(labelsize=23)
        plt.axvline(x=xtaldeg,color="red",linestyle="-",lw=5)
        plt.axhline(y=1.8,color="blue", linestyle="--")
        plt.axhline(y=5.0,color="red", linestyle="--")
    axs[rowidx][3].axvline(x=bestnofilter_deg,color="green",linestyle="-",lw=5)
    axs[rowidx][4].axvline(x=bestfilter_deg,color="green",linestyle="-",lw=5)
    axs[rowidx][3].axvline(x=topnofilter_deg,color="purple",linestyle="-",lw=5)
    axs[rowidx][4].axvline(x=topfilter_deg,color="purple",linestyle="-",lw=5)
    axs[rowidx][1].set_ylim(0,100)
    axs[rowidx][1].set_title(cate,font1)
    axs[rowidx][2].set_ylim(0,100)
    axs[rowidx][2].set_title(f"{len(degstatistics[1])} / {len(degstatistics[0])}" + ", " + format(len(degstatistics[1])/len(degstatistics[0]), '.2%'),font1)
    axs[rowidx][3].set_title(f"nofilter",font1)
    axs[rowidx][4].set_title(f"filter, {len(degstatistics[3])} / {len(degstatistics[2])}" + ", " + format(len(degstatistics[3])/len(degstatistics[2]), '.2%'),font1)

def dealdata(qmgopath, inpdbid, i, xtalligsdf, sdffiles, bestnofiltersdf, bestfiltersdf, topnofiltersdf, topfiltersdf, TS, TQ,xmlfile):
    optlogpath = Path(f"{qmgopath}/{inpdbid}_{i}")
    smarts = TS[i]
    torquartet = TQ[i]
    bincounts = CSD_hist(smarts, xmlfile) # result is the same as that in Torsion Analyzer
    rel_E, func, min_E = QMoptE_list(optlogpath=optlogpath) # QMoptE related infos
    TEU_energies = CSDTEU_list(smarts, xmlfile)

    # mapindex0 = get_map_index(xtalligsdf, sdffiles[0])
    # mapindex1 = get_map_index(xtalligsdf, sdffiles[2])
    mapindex2 = get_map_index(xtalligsdf, bestnofiltersdf)
    # modtorquartet0 = [ mapindex0[i] for i in torquartet ]
    # modtorquartet1 = [ mapindex1[i] for i in torquartet ]
    modtorquartet2 = [ mapindex2[i] for i in torquartet ]
    print(f"sel @/serialNumber={torquartet[0]+1} @/serialNumber={torquartet[1]+1} @/serialNumber={torquartet[2]+1} @/serialNumber={torquartet[3]+1}")
    # print(f"sel @/serialNumber={modtorquartet0[0]+1} @/serialNumber={modtorquartet0[1]+1} @/serialNumber={modtorquartet0[2]+1} @/serialNumber={modtorquartet0[3]+1}")
    # print(f"sel @/serialNumber={modtorquartet1[0]+1} @/serialNumber={modtorquartet1[1]+1} @/serialNumber={modtorquartet1[2]+1} @/serialNumber={modtorquartet1[3]+1}")
    print(f"sel @/serialNumber={modtorquartet2[0]+1} @/serialNumber={modtorquartet2[1]+1} @/serialNumber={modtorquartet2[2]+1} @/serialNumber={modtorquartet2[3]+1}")
    # print()

    sdffunc0 = partial(sdf_ang_list, torquartet=torquartet)
    # sdffunc1 = partial(sdf_ang_list, torquartet=modtorquartet0)
    # sdffunc2 = partial(sdf_ang_list, torquartet=modtorquartet1)
    sdffunc3 = partial(sdf_ang_list, torquartet=modtorquartet2)
    
    # degstatistics = list(map(sdffunc1, sdffiles[0:2])) + list(map(sdffunc2, sdffiles[2:]))
    
    xtaldeg = sdffunc0(xtalligsdf)[0]
    bestnofilter_deg = sdffunc3(bestnofiltersdf)[0]
    bestfilter_deg = sdffunc3(bestfiltersdf)[0]
    topnofilter_deg = sdffunc3(topnofiltersdf)[0]
    topfilter_deg = sdffunc3(topfiltersdf)[0]

    degstatistics=[]
    refmol = Chem.SDMolSupplier(xtalligsdf)[0]
    canonical_mol1 = Chem.MolFromSmiles(Chem.MolToSmiles(refmol,isomericSmiles=False))
    for sdffile in sdffiles:
        degstatistic = []
        mols = Chem.SDMolSupplier(sdffile)
        for i, mol in enumerate(mols):
            canonical_mol2 = Chem.MolFromSmiles(Chem.MolToSmiles(mol,isomericSmiles=False))
            mapindex_ = dict(zip(
                list(map(int, refmol.GetProp("_smilesAtomOutputOrder")[1:-2].split(","))),
                list(map(int, mol.GetProp("_smilesAtomOutputOrder")[1:-2].split(","))),
                ))
            modquartet = [ mapindex_[i] for i in torquartet ]
            idxr1,idxr2,idxr3,idxr4 = modquartet
            deg = rdMolTransforms.GetDihedralDeg(mol.GetConformer(),idxr1,idxr2,idxr3,idxr4)
            degstatistic.append(deg)
        degstatistics.append(degstatistic)
    print(len(degstatistics[0]),len(degstatistics[1]),len(degstatistics[2]),len(degstatistics[3]))

    return smarts,bincounts,func,min_E,rel_E,TEU_energies,degstatistics,xtaldeg,bestnofilter_deg,bestfilter_deg,topnofilter_deg,topfilter_deg

def main(cate):
    print(f">>> Dealing with {cate}")
    inpdbid, samptype, maxconf, readtype = cate.split("-")
    # csvfile = "/pubhome/qcxia02/work/confgen/src/rmsd/crystalligsdf/xtal-lig.uni_Torsion_Strain.csv"
    csvfile = "./xtal-lig.uni_Torsion_Strain.csv" if not readtype == "sani" else "./xtal-lig.uni_sani_Torsion_Strain.csv"
    df = pd.read_csv(csvfile, header=None)
    # pdbidlist = list(df[0])
    pdbidtorsiondict = maketorsiondict(csvfile) # dict[pdbid] = [TS, TQ]
    imgpath = Path("/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/strain/strainimgs") / inpdbid
    TS,TQ = pdbidtorsiondict[inpdbid]

    if samptype != "TLDR":
        cmpdpath = f"/pubhome/qcxia02/work/confgen/compounds/coreset/{maxconf}/{inpdbid}"
        os.system(f"unicon -i {cmpdpath}/conformer.{inpdbid}-{samptype}.fixed.mol2 -o /tmp/{inpdbid}-{samptype}-{maxconf}-{readtype}.sdf &> /dev/null")
        os.system(f"unicon -i {cmpdpath}/conformer.{inpdbid}-{samptype}.fixed.38sani.mol2 -o /tmp/{inpdbid}-{samptype}-{maxconf}-38sani-{readtype}.sdf &> /dev/null")
        dockpath = f"/pubhome/qcxia02/work/confgen/dock/{inpdbid}/blastermaster"
        if not imgpath.exists():
            imgpath.mkdir()
        rownum = len(TS)
        # fig, axs = plt.subplots(rownum, 4, figsize=(60, 10*rownum), sharex=True, sharey=True)
        fig, axs = plt.subplots(rownum, 5, figsize=(60, 10*rownum), sharex=True)
        for i in range(len(TS)):
            xtalligsdf = f"{xtalligpath}/{inpdbid}.sdf"
            sdffiles = [f"/tmp/{inpdbid}-{samptype}-{maxconf}-{readtype}.sdf", f"/tmp/{inpdbid}-{samptype}-{maxconf}-38sani-{readtype}.sdf", f"{dockpath}/docking{samptype}_c{maxconf}b10/test.sdf",f"{dockpath}/docking{samptype}_c{maxconf}b10d38{readtype}/test.sdf"]
            bestnofiltersdf = f"{rmsdsdfpath}/{inpdbid}.{samptype}_c{maxconf}b10.BestRMSDPose.sdf"
            bestfiltersdf = f"{rmsdsdfpath}/{inpdbid}.{samptype}_c{maxconf}b10d38{readtype}.BestRMSDPose.sdf"
            topnofiltersdf = f"{rmsdsdfpath}/{inpdbid}.{samptype}_c{maxconf}b10.TopScorePose.sdf"
            topfiltersdf = f"{rmsdsdfpath}/{inpdbid}.{samptype}_c{maxconf}b10d38{readtype}.TopScorePose.sdf"
            try:
                smarts,bincounts,func,min_E,rel_E,TEU_energies,degstatistics,\
                    xtaldeg,bestnofilter_deg,bestfilter_deg,topnofilter_deg,topfilter_deg = \
                        dealdata(
                            qmgopath, inpdbid, i, xtalligsdf, sdffiles, bestnofiltersdf, \
                                bestfiltersdf, topnofiltersdf, topfiltersdf, TS, TQ,xmlfile
                                )
                # plotall(
                    # smarts,bincounts,func,min_E,rel_E,TEU_energies,degstatistics,xtaldeg,cate=f"{inpdbid}-{samptype}_c{maxconf}b10-{i}",imgpath=imgpath
                # )
                plotallall(
                    smarts,bincounts,func,min_E,rel_E,TEU_energies,degstatistics,xtaldeg,\
                        bestnofilter_deg,bestfilter_deg,topnofilter_deg,topfilter_deg,\
                            cate=f"{inpdbid}-{samptype}_c{maxconf}b10-{i}",\
                                imgpath=imgpath,rowidx=i,axs=axs
                )
            except AttributeError: # Probe torsion smarts does not have a histogram_converted (approximate)
                print(f">>> No 'histogram_converted' was found in Torsion Library v2.1 of {TS[i]}")
                print(f">>> Plot of the corresponding row is left as null.")
            except OSError:
                print(f">>> test.sdf file is null")

        plt.savefig(imgpath / f"{inpdbid}-{samptype}_c{maxconf}b10{readtype}.png",facecolor = "white", transparent = False)
        plt.show()
        plt.close()
        os.system(f"rm /tmp/{inpdbid}-{samptype}-{maxconf}-{readtype}.sdf")
    else:
        cmpdpath = f"/pubhome/qcxia02/work/confgen/compounds/coreset/tldr/valid/{inpdbid}"
        os.system(f"unicon -i {cmpdpath}/{inpdbid}.mol2 -o /tmp/{inpdbid}-{samptype}-{readtype}.sdf &> /dev/null")
        os.system(f"unicon -i {cmpdpath}/{inpdbid}.38sani.mol2 -o /tmp/{inpdbid}-{samptype}-38sani-{readtype}.sdf &> /dev/null")

        dockpath = f"/pubhome/qcxia02/work/confgen/dock/{inpdbid}/blastermaster"
        if not imgpath.exists():
            imgpath.mkdir()
        rownum = len(TS)
        fig, axs = plt.subplots(rownum, 5, figsize=(60, 10*rownum), sharex=True)

        for i in range(len(TS)):
            xtalligsdf = f"{xtalligpath}/{inpdbid}.sdf"
            sdffiles = [f"/tmp/{inpdbid}-{samptype}-{readtype}.sdf", f"/tmp/{inpdbid}-{samptype}-38sani-{readtype}.sdf", f"{dockpath}/docking{samptype}b10/test.sdf",f"{dockpath}/docking{samptype}b10d38{readtype}/test.sdf"] # test.sdf obabeled from test.mol2.gz
            bestnofiltersdf = f"{rmsdsdfpath}/{inpdbid}.{samptype}b10.BestRMSDPose.sdf"
            bestfiltersdf = f"{rmsdsdfpath}/{inpdbid}.{samptype}b10d38{readtype}.BestRMSDPose.sdf"
            topnofiltersdf = f"{rmsdsdfpath}/{inpdbid}.{samptype}b10.TopScorePose.sdf"
            topfiltersdf = f"{rmsdsdfpath}/{inpdbid}.{samptype}b10d38{readtype}.TopScorePose.sdf"

            try:
                smarts,bincounts,func,min_E,rel_E,TEU_energies,degstatistics,\
                    xtaldeg,bestnofilter_deg,bestfilter_deg,topnofilter_deg,topfilter_deg = \
                        dealdata(
                            qmgopath, inpdbid, i, xtalligsdf, sdffiles, bestnofiltersdf, \
                                bestfiltersdf, topnofiltersdf, topfiltersdf, TS, TQ, xmlfile)
                # plotall(
                    # smarts,bincounts,func,min_E,rel_E,TEU_energies,degstatistics,xtaldeg,cate=f"{inpdbid}-{samptype}b10-{i}",imgpath=imgpath
                # )
                plotallall(
                    smarts,bincounts,func,min_E,rel_E,TEU_energies,degstatistics,xtaldeg,\
                        bestnofilter_deg,bestfilter_deg,topnofilter_deg,topfilter_deg,\
                            cate=f"{inpdbid}-{samptype}b10-{i}",\
                                imgpath=imgpath,rowidx=i,axs=axs
                )
            except AttributeError: # Probe torsion smarts does not have a histogram_converted (approximate)
                print(f">>> No 'histogram_converted' was found in Torsion Library v2.1 of {TS[i]}")
                print(f">>> Plot of the corresponding row is left as null.")
            except OSError:
                print(f">>> test.sdf file is null")

        plt.savefig(imgpath / f"{inpdbid}-{samptype}b10{readtype}.png",facecolor = "white", transparent = False)
        plt.show()
        plt.close()

        os.system(f"rm /tmp/{inpdbid}-{samptype}-{readtype}.sdf")
    print(f">>> Finished with {cate}")

if __name__ == "__main__":
    # xmlfile = "/pubhome/qcxia02/git-repo/db2_converter/strain/TL_3.0_VERSION_6.xml"
    xmlfile = "/pubhome/qcxia02/git-repo/db2_converter/strain/TL_2.1_VERSION_6.xml"
    xtalligpath = "/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/pdbbind_mol2s"
    rmsdsdfpath = "/pubhome/qcxia02/work/confgen/src/rmsd/result_rmsd"
    qmgopath = "/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/outputs/PDB/4QMGO"

    cates = []
    # for inpdbid in pdbidlist:
    # for inpdbid in ["1nc3","2cet","2v00","2wtv","3fur","3gc5","4w9c","5dwr"]:
    # for inpdbid in ["2cet","2v00","2wtv","3fur","3gc5","4w9c","5dwr"]:
    # for inpdbid in ["2v00","2wtv","3fur","3gc5","4w9c","5dwr"]:
    for inpdbid in ["2r9w","4w9i","4ty7","4gkm","4jia","4e6q","2vvn","3b27","2w4x","3ryj","4jxs","4bkt","2xii","2xj7","3dd0","2fvd","2w66","2cet","2v00","2wtv","3fur","3gc5","4w9c","5dwr"]:
    # for inpdbid in ["3fur"]:
        for readtype in ["", "sani", "qm"]:
        # for readtype in ["", "sani"]:
        # for readtype in ["sani"]:
        # for readtype in ["qm"]:
            # for samptype in ["C"]:
            for samptype in ["C", "B", "R", "T"]:
                for maxconf in [1000]:
                # for maxconf in [100,250,1000]:
                    cates.append(f"{inpdbid}-{samptype}-{maxconf}-{readtype}")
            cates += [f"{inpdbid}-TLDR--{readtype}"]
        
    # for cate in cates:
        # main(cate)
    with Pool(20) as pool: # cpu cores
        pool.map(main,cates)