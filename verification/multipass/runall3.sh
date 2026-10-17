#!/bin/bash
cd /home/claude/repro
# wait for the in-flight EGFR seed-0 run to finish, then take over with a cheaper order


# fast targets first, all 5 seeds
for s in 0 1 2 3 4; do
  for t in "NT5C2 25 4" "SMO 31 7" "ALK 33 4" "KIT 23 9"; do
    d=fold_$(echo $t|tr ' ' '_')_s$s
    [ -f $d/scores.csv ] || python3 folds.py $t $s >> runall3.log 2>&1
  done
  echo "FAST SEED $s DONE" >> runall3.log
done
# EGFR full: 2 seeds only (seed 0 already running/done)
for s in 0 1; do
  [ -f fold_EGFR_37_30_s$s/scores.csv ] || python3 folds.py EGFR 37 30 $s >> runall3.log 2>&1
done
echo "ALL DONE" >> runall3.log
