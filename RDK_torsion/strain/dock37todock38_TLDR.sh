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

for samptype in "TLDR";do
    for readtype in "qm"; do
        indb2gz=$(ls /pubhome/qcxia02/work/confgen/compounds/coreset/tldr/valid/$inpdbid/$inpdbid.*.db2.gz | tail -1) # dock37 db2
        indb2path=/pubhome/qcxia02/work/confgen/compounds/coreset/tldr/valid/$inpdbid
        indb2gzfile=$(basename $indb2gz)
        tmpdb2gz=/pubhome/qcxia02/work/confgen/compounds/coreset/tldr/valid/$inpdbid/tmp.db2.gz
        cp /pubhome/qcxia02/work/confgen/compounds/coreset/tldr/allbuild3d/$indb2gzfile $tmpdb2gz
        outdb2gz=$indb2path/$(basename $indb2gz .db2.gz).38$readtype.db2.gz # dock38 db2
        rootpath=/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/outputs/PDB
        xtalsdfpath=/pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/rdkit_Pfrag/pdbbind_mol2s

        gunzip $tmpdb2gz
        indb2file=$indb2path/$(basename $tmpdb2gz .gz)
        /pubhome/qcxia02/db22mol2/db2_to_mol2_teb.py $indb2file $indb2path/$inpdbid-$samptype
        gzip $indb2file
        
        for i in $(ls $indb2path/$inpdbid-$samptype.*.mol2 | sort -k2n -t .); do
            cat $i >> $indb2path/$inpdbid-$samptype.mol2
        done

        obabel -imol2 $indb2path/$inpdbid-$samptype.mol2 -osdf -O $indb2path/$inpdbid-$samptype.db2.sdf
        rm $indb2path/$inpdbid-$samptype.*mol2

        python /pubhome/qcxia02/git-repo/TorsionNet/RDK_torsion/strain/calc_strain.py \
        --sdf $indb2path/$inpdbid-$samptype.db2.sdf \
        --refsdf $xtalsdfpath/$inpdbid.sdf \
        --rootpath $rootpath \
        --outpath $indb2path

        csvfile=$indb2path/summary.csv
        /pubhome/qcxia02/db22mol2/dock37todock38.py $indb2path/$(basename $tmpdb2gz) $outdb2gz $csvfile
        rm $csvfile $indb2path/$inpdbid-$samptype.db2.sdf $indb2path/$(basename $tmpdb2gz)
    # break
    done
# break
done
