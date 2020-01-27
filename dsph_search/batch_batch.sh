#!/usr/bin/bash

rm *.slurm
mkdir region_list
python3 the_search.utils.generate_full_sky_cones(3.16)
python3 slurm_maker_region.py

for file in *.slurm;
do
  sbatch ${file}
done

rm *.slurm
