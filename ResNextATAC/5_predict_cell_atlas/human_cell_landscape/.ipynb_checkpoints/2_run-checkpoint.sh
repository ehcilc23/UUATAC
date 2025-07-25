#!/bin/bash 
#SBATCH --job-name=pred_cell_atlas
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --output=%J.out
#SBATCH --error=%J.err

python 2_impute.py --region Region.pred.143w.Human.input.txt

