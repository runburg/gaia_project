#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Make slurm files.

Author: Jack Runburg
Date: 27-09-2019 13:10


"""

coords_list = [(260, 60), (100, -50), (0, -60), (341, 5.7), (37, -80), (159, 51), (152, 10), (15, -35)]
region_radius = 5
# 15 degree close to galactic plane takes ~60 min
num_cones = 10000000
radii = [1.5, 1.0, 0.5]

for num, coords in enumerate(coords_list):
    args = [*coords, region_radius, num_cones, *radii]
    args_string = " ".join([str(arg) for arg in args])

    with open(f'./searchslurm{num}.slurm', 'w') as outfile:
        outfile.write(f'''#!/bin/bash
#SBATCH --job-name=the_search_{num}
#SBATCH --partition=shared
## 3 day max run time for public partitions, except 4 hour max runtime for the sandbox partition
#SBATCH --time=02-00:00:00 ## time format is DD-HH:MM:SS
## task-per-node x cpus-per-task should not typically exceed core count on an individual node
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G ## max amount of memory per node you require

##SBATCH --core-spec=0 ## Uncomment to allow jobs to request all cores on a node

#SBATCH --error=error_{num}_%A.err ## %A - filled with jobid
#SBATCH --output=out_{num}_%A.out ## %A - filled with jobid
## Useful for remote notification
##SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE,TIME_LIMIT_80
##SBATCH --mail-user=runburg@hawaii.edu

source ~/.bash_profile

module load lang/Python/3.7.2-intel-2018.5.274

pip3 install --user numpy
pip3 install --user matplotlib
pip3 install --user astropy
pip3 install --user astroquery
echo "Finished installing packages"
''')
        outfile.write(f"python3 -vu main.py {args_string}")
