import sys
from pathlib import Path
from ase.units import Hartree, kcal, mol

from GOutils import sdf2xtbGOinp_sdf, sdf2SPorcainp, sdf2GOorcainp,SPout2energy,GOout2energy, SPorcainp2GOorcainp

rootpath = Path(sys.argv[1])
startpath = rootpath / "0start"
MMscanpath = rootpath / "1MMscans"
xtbGOpath = rootpath / "2xtbGO"
QMSPpath = rootpath / "3QMSP"
QMGOpath = rootpath / "4QMGO"

subnames = [ item.name for item in QMGOpath.iterdir()]

for subname in subnames:
    outpath = QMGOpath / subname
    outfiles = outpath.glob("*.log")
    for outfile in outfiles:
        with open(QMGOpath / "QMGO.csv",'a') as f:
            energy = GOout2energy(outfile)
            f.write(f"{outfile.name[:-4] + '.orcainp.xyz'},{energy*Hartree/(kcal/mol)}\n")
