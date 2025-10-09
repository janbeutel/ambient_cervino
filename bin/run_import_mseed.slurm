#!/bin/bash
#SBATCH --job-name=cervino_import_mseed
#SBATCH --output=/scratch/jan.beutel/slurm_logs/cervino_import_mseed_%A_%a.log
#SBATCH --error=/scratch/jan.beutel/slurm_logs/cervino_import_mseed_%A_%a.err
##SBATCH --array=0-32%6         # 11 years × 3 channels = 33 tasks
#SBATCH --array=0-23%4          # 8 years × 3 channels = 24 tasks
###SBATCH --time=12:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=2
#SBATCH --partition=IFIall

# Define ranges
#YEARS=(2015 2016 2017 2018 2019 2020 2021 2022 2023 2024 2025)
#STATION=MH36
YEARS=(2015 2016 2017 2018 2019)
STATION=MH38
#YEARS=(2018 2019 2020 2021 2022 2023 2024 2025)
#STATION=MH44
#YEARS=(2019 2020 2021 2022 2023 2024 2025)
#STATION=MH48
#YEARS=(2019 2020)
#STATION=MH54
CHANNELS=(EHE.D EHN.D EHZ.D)

# Map task ID to year/channel
YEAR_IDX=$(( SLURM_ARRAY_TASK_ID / ${#CHANNELS[@]} ))
CHAN_IDX=$(( SLURM_ARRAY_TASK_ID % ${#CHANNELS[@]} ))

YEAR=${YEARS[$YEAR_IDX]}
CHANNEL=${CHANNELS[$CHAN_IDX]}

echo "SLURM running import for YEAR=$YEAR, STATION=$STATION, CHANNEL=$CHANNEL"

# Run MATLAB with year + channel arguments
# matlab -batch "run_import_mseed($YEAR, $STATION, '$CHANNEL')"
matlab -batch "CERVINO_import_mseed_orari($YEAR, '$STATION', '$CHANNEL')"

echo "done."