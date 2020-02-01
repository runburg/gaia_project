#!/usr/bin/bash

rm *.slurm
mkdir region_candidates

echo $"# previous run\n" >> successful_candidates_previous.txt
cat successful_candidates.txt >> successful_candidates_previous.txt
echo "# regions with potential dwarf candidates" > successful_candidates.txt

rm -rf region_list/
mkdir region_list

python3 slurm_maker_region.py

for file in *.slurm;
do
  sbatch ${file}
done

rm *.slurm
