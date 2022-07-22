#$ -S /bin/bash
#$ -q benz
#$ -pe benz 1
#$ -V
#$ -o log
#$ -e log
#$ -j y
#$ -cwd

inpdbid=$inpdbid

conda activate basic
strainpath="/pubhome/qcxia02/git-repo/db2_converter/strain"
tmppath="/tmp/qcxia02/strain"
test -d $tmppath || mkdir -p $tmppath

for samptype in C B R T; do # Conformator, BCL and Rdkit
# for samptype in C; do # Conformator, BCL and Rdkit
    for max_conf in 100 250 1000; do
        # for readtype in "sani" ""; do
        for readtype in "qm"; do
            cate=${inpdbid}-${samptype}
            indb2gz=/pubhome/qcxia02/work/confgen/compounds/coreset/$max_conf/$inpdbid/${cate}.db2.gz # dock37 db2
            outdb2gz=/pubhome/qcxia02/work/confgen/compounds/coreset/$max_conf/$inpdbid/${cate}.38$readtype.db2.gz # dock38 db2
            rootpath=/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/outputs/PDB
            xtalsdfpath=/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/pdbbind_mol2s
            inname=$cate
        
            python /pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/strain/calc_strain.py \
            --sdf \
            --refsdf $xtalsdfpath/$inpdbid.sdf \
            --rootpath $rootpath
            csvfile=$rootpath/summary.csv
            /pubhome/qcxia02/db22mol2/dock37todock38.py $indb2gz $outdb2gz $csvfile
            rm $csvfile $tmppath/${cate}.db2.gz
        # break
        done
    # break
    done
# break
done
