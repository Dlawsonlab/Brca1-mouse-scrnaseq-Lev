#!/bin/bash
#SBATCH -A dalawson_lab
#SBATCH --job-name cellplex
#SBATCH --nodes=1
#SBATCH --mem=64gb
#SBATCH --partition=standard

# load required module

module load cellranger/7.0.1

RESULT_PATH="/share/crsp/lab/dalawson/tzhuravl/CellPlex_Tumors/Cellranger"

cd $RESULT_PATH

echo "started"

cellranger multi --id=nMFP_plus_pMFP_cellplex5 --csv=/share/crsp/lab/dalawson/tzhuravl/CellPlex_Tumors/config_csv/cellplex-normal-premalig-042823.csv

echo "finished"
