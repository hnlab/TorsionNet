import os
import sys

from pathlib import Path


OBABELEXE = "/usr/bin/obabel"

def sdf2xtbGOinp_xyz(sdffile, torsion_quartet, outpath="xtbGO"):
    outpath = Path(outpath)
    if not outpath.exists():
        outpath.mkdir()

    # trans sdf to xyz
    sdffile = Path(sdffile)
    xyzfile = str(sdffile.name)+'.xyz'
    os.system(f"{OBABELEXE} -isdf {str(sdffile)} -oxyz -O {str(outpath / xyzfile)}")
    if (outpath / xyzfile).exists():
        with open(outpath / xyzfile, 'a') as f:
            f.write(f"""
$constrain
force constant=10.0
$end

$constrain
dihedral: {torsion_quartet}, auto
$end
""")
    else:
        print("no file found")

    return  (outpath / xyzfile)


def xyz2SPorcainp(xyzfile, taskline, chg, mult, outpath=""):
    """
    ::input: sdffile (should be absolute path)
    """
    orcainpfile = f"{xyzfile}.orcainp" if not outpath else str(Path(outpath) / f"{Path(xyzfile).name}.orcainp")
    os.system(f"{OBABELEXE} -ixyz {xyzfile} -oorcainp -O {orcainpfile}")
    lines = Path(orcainpfile).read_text().split("\n")
    lines[2] = taskline
    lines[3] = f"* xyz {chg} {mult}"
    Path(orcainpfile).write_text('\n'.join(lines))

    return Path(orcainpfile)

def xyz2GOorcainp(xyzfile, taskline, chg, mult, constrained=False, outpath="", *args):
    """
    ::input: xyzfile (should be absolute path)
    """
    orcainpfile = f"{xyzfile}.orcainp" if not outpath else str(Path(outpath) / f"{Path(xyzfile).name}.orcainp")
    os.system(f"{OBABELEXE} -ixyz {xyzfile} -oorcainp -O {orcainpfile}")
    lines = Path(orcainpfile).read_text().split("\n")
    lines[2] = taskline
    lines[3] = f"* xyz {chg} {mult}"
    if constrained:
        if args:
            torsion_quartet = args[0]
        print("Constrained Optimization is Turned On.")
        orca_torsionatoms = torsion_quartet
        constrain_lines = "%geom\nConstraints\n"+"{D " + orca_torsionatoms +" C}\n" + "end\n" + "end\n"
        lines = "\n".join(lines[:3]) + "\n" + constrain_lines + "\n".join(lines[3:])
    else:
        print("Constrained Optimization is Turned Off.")

    Path(orcainpfile).write_text(lines)
    return Path(orcainpfile)


def SPout2energy(orcaoutfile):
    """
    Specifically for SP out file
    """
    energy=float(0) # 0 is high, which means we skip this conformation

    lines = Path(orcaoutfile).read_text().split("\n")
    for i,line in enumerate(lines):
        if "FINAL SINGLE POINT ENERGY" in line:
            energy = float(line.strip().split()[-1])
            break 
            
    return energy

def GOout2energy(orcaoutfile):
    """
    Specifically for GeomOpt out file
    """
    lines = Path(orcaoutfile).read_text().split("\n")
    for i,line in enumerate(lines):
        # if re.search("*** OPTIMIZATION RUN DONE ***", line):
        if "*** OPTIMIZATION RUN DONE ***" in line:
            energy = float(lines[i-3].split()[4])
    return energy


# if True:
#     sdffile = "/pubhome/qcxia02/cpds.nonbinders200.0_1.sdf"
#     torsion_quartet = "4,6,7,8"
#     outpath = "tmp"

#     sdf2xtbGOinp(sdffile, torsion_quartet, outpath)