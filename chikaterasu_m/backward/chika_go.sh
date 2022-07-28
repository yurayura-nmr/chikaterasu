mkdir -p final_results

./aru_initram.sh -f 0p55_nm_cyclic-CG_to_backcalc.gro -o aa_charmm.gro -to charmm36 -p ./topol.top 
mv aa_charmm.gro final_results/0p55_nm_backcalculated.gro

./aru_initram.sh -f 0p90_nm_cyclic-CG_to_backcalc.gro -o aa_charmm.gro -to charmm36 -p ./topol.top 
mv aa_charmm.gro final_results/0p90_nm_backcalculated.gro

./aru_initram.sh -f 1p30_nm_cyclic-CG_to_backcalc.gro -o aa_charmm.gro -to charmm36 -p ./topol.top 
mv aa_charmm.gro final_results/1p30_nm_backcalculated.gro

./aru_initram.sh -f 0p55_nm_noncyclic-CG_to_backcalc.gro -o aa_charmm.gro -to charmm36 -p ./topol.top 
mv aa_charmm.gro final_results/0p55_nm_backcalculated_noncyclic.gro

./aru_initram.sh -f 0p90_nm_noncyclic-CG_to_backcalc.gro -o aa_charmm.gro -to charmm36 -p ./topol.top 
mv aa_charmm.gro final_results/0p90_nm_backcalculated_noncyclic.gro

./aru_initram.sh -f 1p30_nm_noncyclic-CG_to_backcalc.gro -o aa_charmm.gro -to charmm36 -p ./topol.top 
mv aa_charmm.gro final_results/1p30_nm_backcalculated_noncyclic.gro

./aru_initram.sh -f 2p35_nm_noncyclic-CG_to_backcalc.gro -o aa_charmm.gro -to charmm36 -p ./topol.top 
mv aa_charmm.gro final_results/2p35_nm_backcalculated_noncyclic.gro

scp -r final_results mayuyu:Dropbox/
