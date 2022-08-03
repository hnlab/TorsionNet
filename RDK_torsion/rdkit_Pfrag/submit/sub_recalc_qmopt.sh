conda activate basic
source /pubhome/qcxia02/work/0benchmarks/source_orca.sh
# '''
# qsub_anywhere.py -c "source sub_recalc_qmopt.sh ZINC000170811339_7 45" -j tmp -N ZINC000170811339_7_45 -n 8 -q benz --qsub_now
# '''
# rootpath=/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/outputs/ShuoGu_AmpC_44
# rootpath=/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/outputs/dude_AmpC_259
rootpath=/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/outputs/ShuoGu_D4_256


subdir=$1
angle=$2

python /pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/recalc_qmopt.py \
--rootpath $rootpath \
--subdir $subdir \
--angle $angle \
&> /pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/logs/$subdir-$angle-re-qmopt.log
