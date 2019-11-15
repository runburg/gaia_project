#!/usr/bin/bash

rm *.slurm
python3 slurm_maker_region.py

for file in *.slurm;
do
  sbatch ${file}
done

<<<<<<< HEAD
# rm *.slurm
=======
rm *.slurm
>>>>>>> 45cf07ad8e4e00f3e161b16e169b2f73c20415b3
