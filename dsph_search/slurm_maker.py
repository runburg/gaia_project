#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Make slurm files.

Author: Jack Runburg
Date: 27-09-2019 13:10


"""

num_cones = 1000000
num_per_file = 4000
plot = False

for num in range(num_cones//num_per_file):
    point_start = num * num_per_file
    point_end = (num+1) * num_per_file
    if point_end > num_cones:
        point_end = None

    with open(f'./searchslurm{num}.slurm', 'w') as outfile:
        outfile.write(f'''#!/bin/bash
#SBATCH --job-name=the_search_{num}
#SBATCH --partition=shared
## 3 day max run time for public partitions, except 4 hour max runtime for the sandbox partition
#SBATCH --time=01-00:00:00 ## time format is DD-HH:MM:SS
## task-per-node x cpus-per-task should not typically exceed core count on an individual node
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6400 ## max amount of memory per node you require

##SBATCH --core-spec=0 ## Uncomment to allow jobs to request all cores on a node

#SBATCH --error=%A.err ## %A - filled with jobid
#SBATCH --output=%A.out ## %A - filled with jobid
## Useful for remote notification
##SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE,TIME_LIMIT_80
##SBATCH --mail-user=runburg@hawaii.edu

source ~/.bash_profile

module load lang/Python/3.7.2-intel-2018.5.274

pip3 install --user numpy
pip3 install --user matplotlib
pip3 install --user astropy
pip3 install --user astroquery
pip3 install --user filelock

# cd gaia_project/dsph_search/
python3 -c "import main; main.main(num_cones={num_cones}, point_start={point_start}, point_end={point_end}, plot={plot})"
''')
