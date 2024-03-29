import os
import sys

from pathlib import Path

OBABELEXE = "/usr/bin/obabel"
SDFBOND="{:3d}{:3d}{:3d}{:3d}{:3d}{:3d}{:3d}"
def fmt_single_sdf(sdffile):
    newlines = []
    lines = Path(sdffile).read_text().split("\n")
    atomcount = int(lines[3].split()[0])
    bondcount = int(lines[3].split()[1])
    for i,line in enumerate(lines):
        if i in range(4+atomcount,4+atomcount+bondcount):
            lspl = [int(i) for i in line.split()]
            for _ in range(7 - len(lspl)): # 7 part
                lspl.append(0)
            line = SDFBOND.format(*lspl)
        newlines.append(line)
    return newlines

def sdf2xtbGOinp_sdf(sdffile, torsion_quartet, outpath="xtbGO"):
    outpath = Path(outpath)
    if not outpath.exists():
        outpath.mkdir()

    # trans sdf to xyz
    sdffile = Path(sdffile)
    xtbsdflines = fmt_single_sdf(sdffile)
    newsdffile = outpath / sdffile.name
    with open(newsdffile, "w") as f:
        f.write("\n".join(xtbsdflines))
        f.write("\n\n")
        f.write(f"""
$constrain
force constant=100.0
$end

$constrain
dihedral: {torsion_quartet}, auto
$end
""")
    return newsdffile

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

def sdf2SPorcainp(sdffile, taskline, chg, mult, outpath=""):
    """
    ::input: sdffile (should be absolute path)
    """
    orcainpfile = f"{sdffile}.orcainp" if not outpath else str(Path(outpath) / f"{Path(sdffile).name}.orcainp")
    os.system(f"{OBABELEXE} -isdf {sdffile} -oorcainp -O {orcainpfile}")
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
    maxcycles = None # None means default 3N
    if constrained:
        if len(args) == 1:
            torsion_quartet = args[0]
            
        if len(args) == 2:
            torsion_quartet = args[0]
            maxcycles = args[1]

        print("Constrained Optimization is Turned On.")
        orca_torsionatoms = torsion_quartet
        if maxcycles == None:
            constrain_lines = "%geom\nConstraints\n"+"{D " + orca_torsionatoms +" C}\n" + "end\n" + "end\n"
        else:
            constrain_lines = "%geom\nMaxIter " + str(maxcycles) + "\nConstraints\n"+"{D " + orca_torsionatoms +" C}\n" + "end\n" + "end\n"
        lines = "\n".join(lines[:3]) + "\n" + constrain_lines + "\n".join(lines[3:])
    else:
        print("Constrained Optimization is Turned Off.")

    Path(orcainpfile).write_text(lines)
    return Path(orcainpfile)

def sdf2GOorcainp(sdffile, taskline, chg, mult, constrained=False, outpath="", *args):
    """
    ::input: sdffile (should be absolute path)
    """
    orcainpfile = f"{sdffile}.orcainp" if not outpath else str(Path(outpath) / f"{Path(sdffile).name}.orcainp")
    os.system(f"{OBABELEXE} -isdf {sdffile} -oorcainp -O {orcainpfile}")
    lines = Path(orcainpfile).read_text().split("\n")
    lines[2] = taskline
    lines[3] = f"* xyz {chg} {mult}"
    maxcycles = None # None means default 3N
    if constrained:
        if len(args) == 1:
            torsion_quartet = args[0]
            
        if len(args) == 2:
            torsion_quartet = args[0]
            maxcycles = args[1]

        print("Constrained Optimization is Turned On.")
        orca_torsionatoms = torsion_quartet
        if maxcycles == None:
            constrain_lines = "%geom\nConstraints\n"+"{D " + orca_torsionatoms +" C}\n" + "end\n" + "end\n"
        else:
            constrain_lines = "%geom\nMaxIter " + str(maxcycles) + "\nConstraints\n"+"{D " + orca_torsionatoms +" C}\n" + "end\n" + "end\n"
        lines = "\n".join(lines[:3]) + "\n" + constrain_lines + "\n".join(lines[3:])
    else:
        print("Constrained Optimization is Turned Off.")

    Path(orcainpfile).write_text(lines)
    return Path(orcainpfile)

def SPorcainp2GOorcainp(SPorcainp, taskline, constrained=False, outpath="", *args):
    """
    ::input: SPorcainp (should be absolute path)
    """
    GOorcainpfile = SPorcainp if not outpath else str(Path(outpath) / Path(SPorcainp).name)
    lines = Path(SPorcainp).read_text().split("\n")
    maxcycles = None # None means default 3N
    if constrained:
        if len(args) == 1:
            torsion_quartet = args[0]
            
        if len(args) == 2:
            torsion_quartet = args[0]
            maxcycles = args[1]
    
        print("Constrained Optimization is Turned On.")
        orca_torsionatoms = torsion_quartet
        if maxcycles == None:
            constrain_lines = "%geom\nConstraints\n"+"{D " + orca_torsionatoms +" C}\n" + "end\n" + "end\n"
        else:
            constrain_lines = "%geom\nMaxIter " + str(maxcycles) + "\nConstraints\n"+"{D " + orca_torsionatoms +" C}\n" + "end\n" + "end\n"
    else:
        print("Constrained Optimization is Turned Off.")

    lines = "\n".join(lines[:2]) + "\n" + taskline + constrain_lines + "\n".join(lines[4:])
    Path(GOorcainpfile).write_text(lines)
    return Path(GOorcainpfile)


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
