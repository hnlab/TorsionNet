#!/bin/bash

rootpath=/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag
workpath=/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/submit
outpath=/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/outputs/PDB
# outpath=/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/outputs/tl-tn-r2scan
MMscanpath=$outpath/1MMscans

for dir in $(ls $MMscanpath); do
# for dir in $dirs; do
# for dir in `cat /tmp/tmp`; do
name=$(basename $dir)
qsub_anywhere.py -c "sh $workpath/sub_all.sh $outpath $name" -q benz -n 8 -j tmp -N all-$name --qsub_now
done
