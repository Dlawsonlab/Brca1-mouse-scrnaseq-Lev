#!/bin/bash
#SBATCH -A ECOEVO283
#SBATCH --job-name cellplex
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --partition=standard

# load required module

module load cellranger/7.0.1

RESULT_PATH="/share/crsp/lab/dalawson/tzhuravl/CellPlex_Tumors/Cellranger"

cd $RESULT_PATH

echo "started normal"

cellranger multi --id=nMFP_plus_Tumor_599 --csv=/share/crsp/lab/dalawson/tzhuravl/CellPlex_Tumors/config_csv/cellplex-normal-tum-012623.csv

echo "finished normal, starting premalignant"

cellranger multi --id=pMFP_plus_Tumor_633 --csv=/share/crsp/lab/dalawson/tzhuravl/CellPlex_Tumors/config_csv/cellplex-premalig-tum-020723.csv

echo "both are finished"
