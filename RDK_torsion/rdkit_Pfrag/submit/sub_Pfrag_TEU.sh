# mol2path=/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/ffparam
# outpath=/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/outputs/tl-tn-b3lyp
# imgpath=/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/imgs/tl-tn-b3lyp

# mol2path=/pubhome/qcxia02/Downloads/dataset/ShuoGu/D4
# outpath=/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/outputs/ShuoGu_D4_256
# imgpath=/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/imgs/ShuoGu_D4_256

# mol2path=/pubhome/qcxia02/Downloads/dataset/ShuoGu/D4
# outpath=/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/outputs/ShuoGu_D4_256
# imgpath=/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/imgs/ShuoGu_D4_256

# mol2path=/pubhome/qcxia02/git-repo/DL4molcst/scripts/dock_models/3-conf-docks/results/AmpC/compounds/mol2mols
# outpath=/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/outputs/dude_AmpC_259/0start
# imgpath=/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/imgs/dude_AmpC_259

# mol2path=/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/PROTAC_mol2s
# outpath=/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/outputs/PROTAC
# imgpath=/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/imgs/PROTAC

mol2path=/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/pdbbind_mol2s
outpath=/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/outputs/PDB
imgpath=/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/imgs/PDB

srcpath=/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag

conda activate basic

# for file in $(ls $mol2path/*.mol2); do
for file in $(ls $mol2path/*.mol2); do
    mol2file=$(basename $file)
    python $srcpath/rdkit_Pfrag_TEU.py \
    --mol2 $mol2path/$mol2file \
    --outpath $outpath \
    --imgpath $imgpath \
    --rmsd 0.1 \
    --numConfs 20 \
    &>> Pfrag_TEU_log
done
