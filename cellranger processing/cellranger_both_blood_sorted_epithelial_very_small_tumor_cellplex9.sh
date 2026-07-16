#!/bin/bash
#SBATCH -A dalawson_lab
#SBATCH --job-name cellplexranger
#SBATCH --nodes=1
#SBATCH --mem=64gb
#SBATCH --partition=standard

# load required module

module load cellranger/7.0.1

RESULT_PATH="/share/crsp/lab/dalawson/share/1_Lab_sequencing/1_Tatyana/Cellranger"

cd $RESULT_PATH

echo "started"

cellranger multi --id=both_blood_sorted_epithelial_very_small_tumor --csv=/share/crsp/lab/dalawson/share/1_Lab_sequencing/1_Tatyana/config_csv/cellplex_both_blood_sorted_epithelial_very_small_tumor_121923.csv

echo "finished"
