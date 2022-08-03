#!/bin/bash
#$ -S /bin/bash
#$ -V
#$ -notify
##$ -l hostname=!k[121,136].hn.org
#$ -cwd

source /pubhome/qcxia02/work/0benchmarks/source_orca.sh

rootpath=$1
subdir=$2
# python /pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/geom_opt.py \
    # --rootpath $1 \
    # --subdir $2 \
    # --MMopt \
    # &> /pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/logs/$2_xtbopt.log
# 
# python /pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/geom_opt.py \
    # --rootpath $1 \
    # --subdir $2 \
    # --QMSP \
    # &> /pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/logs/$2_qmsp.log

# python /pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/geom_opt.py \
    # --rootpath $1 \
    # --subdir $2 \
    # --QMopt \
    # &> /pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/logs/$2_qmopt.log

python /pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/geom_opt.py \
    --rootpath $rootpath \
    --subdir $subdir \
    --QMopt \
    --method "B3LYP def2-SVP" \
    &> /pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/logs/${subdir}_qmopt.log

