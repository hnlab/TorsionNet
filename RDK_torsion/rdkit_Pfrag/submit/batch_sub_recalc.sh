srcpath=/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/submit

subdir=D4.76_1
angle=150
qsub_anywhere.py -c "source $srcpath/sub_recalc_qmopt.sh $subdir $angle" -q opel -n 8 -j tmp -N $subdir-$angle-re-qmopt --qsub_now

subdir=D4.76_1
angle=-150
qsub_anywhere.py -c "source $srcpath/sub_recalc_qmopt.sh $subdir $angle" -q opel -n 8 -j tmp -N $subdir-$angle-re-qmopt --qsub_now

