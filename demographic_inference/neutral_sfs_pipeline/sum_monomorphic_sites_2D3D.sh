#! /bin/bash
#SBATCH --partition=GPU

cd /data/shared/snigenda/finwhale_projects/fin_genomics/demographic_inference/SFS/Neutral/SFS_projection_Monomorphic

########################
# 2D (pair) sums
########################

out="2D_monomorphic_sums.txt"

echo -e "pair\tmonomorphic_sites" > "$out"

for pair in ENP-GOC ENP-ESP ESP-GOC; do
  awk -v p1=$(echo "$pair" | cut -d- -f1) \
      -v p2=$(echo "$pair" | cut -d- -f2) '
    FNR>1 && $1==p1 && $2==p2 {sum+=$3}
    END {print p1"-"p2, sum+0}
  ' countsOfMonomorphicPassingProjectionThresholds.perPair.{01..21}.txt \
  | tee -a "$out"
done

echo "2D results written to $out"


########################
# 3D (triple) sums
########################

out3D="3D_monomorphic_sums.txt"

echo -e "triple\tmonomorphic_sites" > "$out3D"

for triple in ENP-ESP-GOC; do
  awk -v t1=$(echo "$triple" | cut -d- -f1) \
      -v t2=$(echo "$triple" | cut -d- -f2) \
      -v t3=$(echo "$triple" | cut -d- -f3) '
    FNR>1 && $1==t1 && $2==t2 && $3==t3 {sum+=$4}
    END {print t1"-"t2"-"t3, sum+0}
  ' countsOfMonomorphicPassingProjectionThresholds.perTriple.{01..21}.txt \
  | tee -a "$out3D"
done

echo "3D results written to $out3D"
