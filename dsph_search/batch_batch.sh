#!/usr/bin/bash

rm *.slurm
mkdir region_candidates
python3 slurm_maker_region.py

for file in *.slurm;
do
  sbatch ${file}
done

rm *.slurm
