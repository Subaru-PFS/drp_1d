#!/bin/bash

dataset="batch6"

nodes="
node17.lam.fr
node18.lam.fr
node19.lam.fr
node20.lam.fr
node21.lam.fr
node22.lam.fr
node23.lam.fr
node24.lam.fr
node25.lam.fr
node26.lam.fr
node27.lam.fr
node28.lam.fr
node29.lam.fr
node30.lam.fr
node31.lam.fr
node32.lam.fr
node33.lam.fr
"

for node in $nodes; do
    qsub -q batch -v dataset=$dataset -l host=$node /net/fs/users/rborges/tmp_clean_run.sh
done
