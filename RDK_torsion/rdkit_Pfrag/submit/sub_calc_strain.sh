# mol2path=/pubhome/qcxia02/Downloads/dataset/ShuoGu/AmpC
mol2path=/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/pdbbind_mol2s
# mol2path=/tmp
# mol2path=/pubhome/qcxia02/git-repo/DL4molcst/scripts/dock_models/3-conf-docks/results/AmpC/studies/total/docking/mol2
# rootpath=/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/outputs/ShuoGu_AmpC_44
# rootpath=/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/outputs/dude_AmpC_259
rootpath=/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/outputs/PDB
srcpath=/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag

conda activate basic

for file in $(ls $mol2path/5dwr.mol2); do
# for file in $(cat /tmp/tmp); do
    mol2file=$(basename $file)
    python $srcpath/calc_strain.py --rootpath $rootpath --mol2 $mol2path/$mol2file &>>$rootpath/strain_log
done
