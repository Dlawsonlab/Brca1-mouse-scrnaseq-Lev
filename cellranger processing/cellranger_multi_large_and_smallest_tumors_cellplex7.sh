#!/bin/bash
#SBATCH -A dalawson_lab
#SBATCH --job-name cellplex
#SBATCH --nodes=1
#SBATCH --mem=64gb
#SBATCH --partition=standard

# load required module

module load cellranger/7.0.1

RESULT_PATH="/share/crsp/lab/dalawson/share/1_Lab_sequencing/1_Tatyana/Cellranger"

cd $RESULT_PATH

echo "started"

cellranger multi --id=large_and_smallest_tumors_cellplex7 --csv=/share/crsp/lab/dalawson/share/1_Lab_sequencing/1_Tatyana/config_csv/cellplex_tumors_082223.csv

echo "finished"
