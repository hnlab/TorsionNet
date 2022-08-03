#!/bin/bash
# rootpath=$1
rootpath=/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag
# outpath=/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/outputs/ShuoGu_AmpC_44
# outpath=/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/outputs/ShuoGu_AmpC_44_new
# outpath=/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/outputs/dude_AmpC_259
# outpath=/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/outputs/ShuoGu_D4_256
outpath=/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/outputs/PDB
# outpath=/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/outputs/tl-tn-r2scan
workpath=/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/submit
MMscanpath=$outpath/1MMscans

# MMscanpath=$outpath/1MMscans/tmp
i=1
# for dir in $(ls $MMscanpath); do
for dir in 4w9c_0 4w9c_1 4w9c_2 4w9c_3 4w9c_4 4w9c_5 4w9c_6 4w9c_7 ; do
# for dir in $(cat /tmp/tmp); do
# for dir in "2y6";do
# for dir in $(cat /tmp/sp.txt); do
# for dir in $(cat /tmp/go.txt); do
# for dir in $(cat /tmp/sub_go.txt); do


class=$(expr $i % 4)
name=$(basename $dir)
### XTBGO ###
# test $class = 0 && qsub_anywhere.py -c "source $workpath/sub.sh $outpath $name" -q benz -n 32 -j tmp -N MMopt_2-$name --qsub_now
# test $class = 1 && qsub_anywhere.py -c "source $workpath/sub.sh $outpath $name" -q benz -n 32 -j tmp -N MMopt_2-$name --qsub_now

# test $class = 0 && qsub -v "rootpath=$outpath" -v "subdir=$name" -q benz -pe benz 32 -j y -N xtbgo-$name sub.sh
# test $class = 1 && qsub -v "rootpath=$outpath" -v "subdir=$name" -q benz -pe benz 32 -j y -N xtbgo-$name sub.sh

### QMSP ###
# test $class = 0 && qsub_anywhere.py -c "source $workpath/sub_SP.sh $outpath $name" -q opel -n 8 -j tmp -N QMSP_2-$name --qsub_now
# test $class = 1 && qsub_anywhere.py -c "source $workpath/sub_SP.sh $outpath $name" -q opel -n 8 -j tmp -N QMSP_2-$name --qsub_now
# test $class = 2 && qsub_anywhere.py -c "source $workpath/sub_SP.sh $outpath $name" -q opel -n 8 -j tmp -N QMSP_2-$name --qsub_now
# test $class = 3 && qsub_anywhere.py -c "source $workpath/sub_SP.sh $outpath $name" -q benz -n 8 -j tmp -N QMSP_2-$name --qsub_now

# test $class = 0 && qsub -q opel -l hostname=!k136.hn.org -pe opel 8 -j y -N qmsp-$name sub_SP.sh
# test $class = 1 && qsub -q mazda -l hostname=!k121.hn.org -pe mazda 8 -j y -N qmsp-$name sub_SP.sh

### QMGO ###
qsub_anywhere.py -c "source $workpath/sub_GO.sh $outpath $name" -q benz -n 32 -j tmp -N QMopt-$name --qsub_now

# test $class = 0 && qsub -q opel -l hostname=!k121.hn.org -pe opel 8 -j y -N qmgo-$name sub_GO.sh
# test $class = 1 && qsub -q opel -l hostname=!k136.hn.org -pe mazda 8 -j y -N qmgo-$name sub_GO.sh
# test $class = 0 && qsub -q benz -pe benz 8 -j y -N qmgo-$name -p -100 sub_GO.sh
# test $class = 1 && qsub -q benz -pe benz 8 -j y -N qmgo-$name -p -100 sub_GO.sh

#######################################################################################
### ALL ###
# test $class = 0 && qsub -v "rootpath=$outpath" -v "subdir=$name" -q benz  -pe benz 8 -j y -N all-$name sub_all.sh
# test $class = 1 && qsub -v "rootpath=$outpath" -v "subdir=$name" -q benz  -pe benz 8 -j y -N all-$name sub_all.sh
# 
#######################################################################################

i=$(expr $i + 1)
# fi
done
