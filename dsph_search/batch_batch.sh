#!/usr/bin/bash

rm *.slurm
python3 slurm_maker_region.py

for file in *.slurm;
do
  sbatch ${file}
done

rm *.slurm
