#!/bin/bash
#SBATCH -A dalawson_lab
#SBATCH --job-name cellplex
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=standard

# load required module

module load cellranger/7.0.1

RESULT_PATH="/share/crsp/lab/dalawson/tzhuravl/CellPlex_Metastasis/Cellranger"

cd $RESULT_PATH

echo "started"

cellranger multi --id=lung_liver_ctrl_tumor_cellplex4 --csv=/share/crsp/lab/dalawson/tzhuravl/CellPlex_Metastasis/config_csv/cellplex-lung-liver-mets-ctrl-tumor_041223.csv

echo "finished"
