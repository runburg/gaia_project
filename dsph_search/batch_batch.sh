#!/usr/bin/env bash
python3 slurm_maker.py
rm *.slurm

for file in *.slurm;
do
  sbatch file
done

rm*.slurm
