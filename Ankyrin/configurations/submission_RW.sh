#!/bin/bash
for d in */;
do
cd -- $d
sbatch job_SMD.sh
cd ..
done
