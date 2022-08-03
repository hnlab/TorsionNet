#!/bin/bash
#$ -S /bin/bash
#$ -notify
##$ -l hostname=!k[121,136].hn.org
#$ -cwd
rootpath=$1
subdir=$2

source /pubhome/qcxia02/work/0benchmarks/source_orca.sh
python /pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/geom_opt.py \
    --rootpath $rootpath \
    --subdir $subdir \
    --MMopt \
    &> /pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/logs/${subdir}_xtbopt.log

# python /pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/geom_opt.py \
#     --rootpath $1 \
#     --subdir $2 \
#     --QMSP \
#     &> /pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/logs/$2_qmsp.log
# # 
# python /pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/geom_opt.py \
    # --rootpath $1 \
    # --subdir $2 \
    # --QMopt \
    # &> /pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/logs/$2_qmopt.log
# 
