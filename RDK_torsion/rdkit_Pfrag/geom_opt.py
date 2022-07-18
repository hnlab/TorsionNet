#!/pubhome/qcxia02/miniconda3/envs/basic/bin/python
import os
import sys
import argparse
import pandas as pd
from pathlib import Path, PosixPath

from rdkit import Chem
from ase.units import Hartree, kcal, mol
molunit = mol
sys.path.append("/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/utils")
from GOutils import sdf2xtbGOinp_sdf, sdf2SPorcainp, sdf2GOorcainp,SPout2energy,GOout2energy, SPorcainp2GOorcainp
from utils import getsdfchg

OBABELEXE = "/usr/bin/obabel"
XTBEXE = "/pubhome/qcxia02/miniconda3/envs/basic/bin/xtb"
ORCAEXE = "/pubhome/soft/orca/orca_5_openmpi411/orca"


if __name__ == "__main__":
    parser = argparse.ArgumentParser("geom_opt")
    parser.add_argument("--rootpath", type=PosixPath, help="absolute path for rootpath")
    parser.add_argument("--subdir", type=str, help="subdir name for computation, not absolute path", required=True)
    parser.add_argument("--MMopt", action="store_true", default=False, help="trigger MM optimization")
    parser.add_argument("--QMSP", action="store_true", default=False, help="trigger QM single point")
    parser.add_argument("--QMopt", action="store_true", default=False, help="trigger QM optimization")
    parser.add_argument("--torsionquartet", nargs="+", help="set torsion quartet manually")
    parser.add_argument("--method", type=str, help="QM level, e.g. 'B3LYP def2-SVP",default="r2SCAN-3c")
    args = parser.parse_args()

    # rootpath = Path("/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/outputs")
    rootpath = args.rootpath
    startpath = rootpath / "0start"
    MMscanpath = rootpath / "1MMscans"
    xtbGOpath = rootpath / "2xtbGO"
    QMSPpath = rootpath / "3QMSP"
    QMGOpath = rootpath / "4QMGO"
    
    for path in [MMscanpath,xtbGOpath,QMSPpath,QMGOpath]:
        if not path.exists():
            path.mkdir()

    if args.MMopt:
        # 1 GFN-FF opt
        # MMscanpath = Path(
            # "/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/outpath/1MMscans"
        # ) # rootpath
        # for dir in MMscanpath.iterdir():
            # outpath = Path("/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/outputs/2xtbGO") / dir.name

        subdir = args.subdir
        MMscanpath_sub = MMscanpath /  subdir # subsetpath

        ## get info from sdffile ##
        sdffile = MMscanpath / MMscanpath_sub / (MMscanpath_sub.name + "-1_000.sdf")
        mol = Chem.SDMolSupplier(str(sdffile))[0]
        try:
            torsion_quartet = mol.GetProp("TORSION_ATOMS_FRAGMENT")
        except:
            print("!Notion: using manually given torsion quartet")
            if args.torsionquartet:
                torsion_quartet = " ".join(args.torsionquartet)

        torsion_quartet_add1 = ",".join(
            list([str(int(idx) + 1) for idx in torsion_quartet.split()])
        )
        # totchg = getsdfchg(sdffile)

        # with open(xtbGOpath / "name_chg_tor-quartet.csv", 'a') as f:
            # f.write(f"{MMscanpath_sub.name},{totchg},{torsion_quartet}\n")
        #####################f######

        outpath = xtbGOpath / MMscanpath_sub.name
        if not outpath.exists():
            outpath.mkdir()
        for sdffile in MMscanpath_sub.iterdir():
            mol = Chem.SDMolSupplier(str(sdffile))[0]
            # prepare input file
            inpfile = sdf2xtbGOinp_sdf(sdffile, torsion_quartet_add1, outpath)
            # run xtb opt
            os.system(f"{XTBEXE} {inpfile} --opt --gfnff")
            # xtboptoutfile = "xtblast.xyz" # xtblast.xyz means sth. wrong with xtb opt
            xtboptoutfile = "xtbopt.sdf" # xtblast.xyz means sth. wrong with xtb opt
            xtboptoutfile_rename = sdffile.name + ".opt.sdf"
            # save optimized xyz
            os.system(f"mv {xtboptoutfile} {outpath/xtboptoutfile_rename}")
            os.system(f"rm xtb* gfn* .xtboptok")

        
    # exclude = Path("/pubhome/qcxia02/code-stdin-x6i").read_text().split("\n")
    if args.QMSP:
        # (1.5 r2SCAN-3c SP)
        # os.system("source /pubhome/qcxia02/work/0benchmarks/source_orca.sh")

        if not QMSPpath.exists():
            QMSPpath.mkdir()

        subdir = args.subdir
        # if not subdir in exclude:
        MMscanpath_sub = MMscanpath / subdir  # subsetpath
        xtbGOpath_sub = xtbGOpath / subdir
        sdffile = MMscanpath / MMscanpath_sub / (MMscanpath_sub.name + "-1_000.sdf")
        mol = Chem.SDMolSupplier(str(sdffile))[0]

        try:
            torsion_quartet = mol.GetProp("TORSION_ATOMS_FRAGMENT")
        except:
            print("!Notion: using manually given torsion quartet")
            if args.torsionquartet:
                torsion_quartet = " ".join(args.torsionquartet)

        totchg = getsdfchg(str(sdffile))
        mult = 1

        outpath = QMSPpath / MMscanpath_sub.name
        if not outpath.exists():
            outpath.mkdir()
        # method = "r2SCAN-3c"
        # method = "B3LYP def2-SVP"
        method = args.method

        taskline = f"! {method}\n%pal nprocs 8 end"  # parallel
        for outsdffile in xtbGOpath_sub.iterdir():
            if outsdffile.name.endswith(".opt.sdf"):
                inpfile = sdf2SPorcainp(
                    outsdffile, taskline, chg=totchg, mult=mult, outpath=outpath
                )
                outfile = outpath / (inpfile.name + ".log")
                os.system(f"{ORCAEXE} {inpfile} &> {outfile}")
                os.system(
                    f"rm {outpath / (inpfile.name + '.gbw')} {outpath / (inpfile.name + '_property.txt')} {outpath / (inpfile.name + '.densities')}"
                )

    if args.QMopt:
        # 2 r2SCAN-3c orca opt
        # os.system("source /pubhome/qcxia02/work/0benchmarks/source_orca.sh")

        if not QMGOpath.exists():
            QMGOpath.mkdir()
        subdir = args.subdir

        QMSPpath_sub = QMSPpath / subdir

        MMscanpath_sub = MMscanpath / subdir  # subsetpath
        xtbGOpath_sub = xtbGOpath / subdir
        sdffile = MMscanpath / MMscanpath_sub / (MMscanpath_sub.name + "-1_000.sdf")
        mol_ = Chem.SDMolSupplier(str(sdffile))[0]

        try:
            torsion_quartet = mol_.GetProp("TORSION_ATOMS_FRAGMENT")
        except:
            print("!Notion: using manually given torsion quartet")
            if args.torsionquartet:
                torsion_quartet = " ".join(args.torsionquartet)

        totchg = getsdfchg(str(sdffile))
        mult = 1

        # MMoptfiles = []
        QMSPfiles = []
        # select confs with minimum energy for each angle
        for ang in range(-180, 180, 15):
            files = list(QMSPpath_sub.glob(f'{subdir}-*_{"%03d" % ang}.sdf.opt.sdf.orcainp.log'))
            energies =  list([ SPout2energy(file) for file in files ])
            minfile = files[energies.index(min(energies))]
            # minoptfilename = minfile.name[:-12]
            minspfile = minfile.name[:-4]
            with open(QMGOpath / "QMSP.csv",'a') as f:
                for i in range(len(files)):
                    f.write(f"{files[i].name},{energies[i]*Hartree/(kcal/molunit)}\n")
            
            # MMoptfiles.append(xtbGOpath / subdir / minoptfilename)
            QMSPfiles.append(QMSPpath / subdir / minspfile)


        outpath = QMGOpath / MMscanpath_sub.name
        if not outpath.exists():
            outpath.mkdir()

        # method = "r2SCAN-3c"
        # method = "B3LYP def2-SVP"
        method=args.method
        taskline = f"! {method} opt\n%pal nprocs 8 end"  # parallel

        energies = []
        # for optsdffile in QMSPfiles:
        for qmspfile in QMSPfiles:
            # if optsdffile.name.endswith(".opt.sdf"):
            # inpfile = sdf2GOorcainp(
                # optsdffile,
                # taskline,
                # totchg,
                # mult,
                # True,
                # outpath,
                # torsion_quartet,
            # )

            inpfile = SPorcainp2GOorcainp(
                qmspfile,
                taskline,
                True,
                outpath,
                torsion_quartet,
            )

            outfile = outpath / (inpfile.name + ".log")
            os.system(f"{ORCAEXE} {inpfile} &> {outfile}")
            os.system(
                # f"rm {outpath / (inpfile.name + '.gbw')} {outpath / (inpfile.name + '_property.txt')} {outpath / (inpfile.name + '.densities')} {outpath / (inpfile.name + '.opt')} {outpath / (inpfile.name + '.engrad')} {outpath / (inpfile.name + '_trj.xyz')}"
                f"rm {outpath / (inpfile.name + '.gbw')} {outpath / (inpfile.name + '_property.txt')} {outpath / (inpfile.name + '.densities')} {outpath / (inpfile.name + '.opt')} {outpath / (inpfile.name + '.engrad')}" # not remove trajectory
            )
            try:
                energy = GOout2energy(outfile)
                with open(QMGOpath / "QMGO.csv",'a') as f:
                    # f.write(f"{optsdffile.name + '.orcainp.xyz'},{energy*Hartree/(kcal/molunit)}\n")
                    f.write(f"{qmspfile.name + '.orcainp.xyz'},{energy*Hartree/(kcal/molunit)}\n")

            except UnboundLocalError:
                print(f"Energy cannot be read from orca out file of {qmspfile.name}, \nplease check carefully if there are any problems in orca optimization")
        # 3 Comparison
        # 1) compare different MMscan starting point PES
        # 2) compare MMscan SP-selected and QM opt best
        # 3) compare MMscan SP-selected and single starting point in 1)
