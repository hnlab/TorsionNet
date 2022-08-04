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
            indb2path=$(dirname $(readlink -f $indb2gz))
            outdb2gz=/pubhome/qcxia02/work/confgen/compounds/coreset/$max_conf/$inpdbid/${cate}.38$readtype.db2.gz # dock38 db2
            rootpath=/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/outputs/PDB
            xtalsdfpath=/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/pdbbind_mol2s
            inname=$cate

            gunzip $indb2gz
            indb2file=$indb2path/$(basename $indb2gz .gz)
            /pubhome/qcxia02/db22mol2/db2_to_mol2_teb.py $indb2file $indb2path/$inpdbid-$samptype-$max_conf
            gzip $indb2file
            
            for i in $(ls $indb2path/$inpdbid-$samptype-$max_conf.*.mol2 | sort -k2n -t .); do
                cat $i >> $indb2path/$inpdbid-$samptype-$max_conf.mol2
            done

            obabel -imol2 $indb2path/$inpdbid-$samptype-$max_conf.mol2 -osdf -O $indb2path/$inpdbid-$samptype-$max_conf.db2.sdf
            rm $indb2path/$inpdbid-$samptype-$max_conf.*mol2

            python /pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/strain/calc_strain.py \
            --sdf $indb2path/$inpdbid-$samptype-$max_conf.db2.sdf \
            --refsdf $xtalsdfpath/$inpdbid.sdf \
            --rootpath $rootpath \
            --outpath $indb2path

            csvfile=$indb2path/summary.csv
            /pubhome/qcxia02/db22mol2/dock37todock38.py $indb2gz $outdb2gz $csvfile
            rm $csvfile $indb2path/$inpdbid-$samptype-$max_conf.db2.sdf
        # break
        done
    # break
    done
# break
done
