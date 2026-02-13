# Short testing (1 ps)
./ai2bmd --prot-file ub.pdb --preprocess-dir ub_preprocessed --preeq-steps 0 --sim-steps 1000 --record-per-steps 1

# Production run (10 ns)
./ai2bmd --prot-file ub.pdb --preprocess-dir ub_preprocessed --preeq-steps 0 --sim-steps 10000000 --record-per-steps 1000


# Current structure this folder
# ub.pdb
# ub_preprocessed

# inside ub_preprocessed
# ub-preeq-nowat.pdb  ub-preeq.pdb
