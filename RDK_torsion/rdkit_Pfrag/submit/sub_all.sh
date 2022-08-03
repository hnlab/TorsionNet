#!/bin/bash
#$ -S /bin/bash
#$ -notify
##$ -l hostname=!k[121,136].hn.org
#$ -cwd

source /pubhome/qcxia02/work/0benchmarks/source_orca.sh
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

# python /pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/geom_opt.py \
#     --rootpath $1 \
#     --subdir $2 \
#     --MMopt \
#     --QMSP \
#     --QMopt \
#     --method "B3LYP def2-SVP" \
#     &> /pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/logs/$2_all.log

# r2scan-3c
python /pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/geom_opt.py \
    --rootpath $1 \
    --subdir $2 \
    --MMopt \
    --QMSP \
    --QMopt \
    &> /pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/logs/$2_all.log