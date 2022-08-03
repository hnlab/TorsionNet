# mol2path=/pubhome/qcxia02/Downloads/dataset/ShuoGu/AmpC
srcpath=/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag
# rootpath=/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/outputs/ShuoGu_AmpC_44
# rootpath=/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/outputs/dude_AmpC_259
# rootpath=/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/outputs/ShuoGu_D4_256
rootpath=/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/outputs/PDB


startpath=$rootpath/0start
qmgopath=$rootpath/4QMGO
conda activate basic

# for file in $(ls $startpath); do
# for file in $(ls $qmgopath); do
# for file in $(cat /pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/outputs/ShuoGu_D4_256/head100.txt); do

# for file in `cat $rootpath/go.txt`; do
# for file in `cat /tmp/tmp`; do
for file in $(ls $startpath); do
    subdir=$(basename $file .sdf)
    python $srcpath/write_sdf_func.py --rootpath $rootpath --subdir $subdir &>> $rootpath/func_log
done
