mkdir chika_ai2bmd
cd chika_ai2bmd

cp from_last_step.pdb ./protein-pre-eq.pdb
../aibmd --prot-file protein-pre-eq.pdb --sim-steps 10000 --gpus all --temp-k 292 --preeq-steps 0
