#!/bin/bash
#SBATCH --job-name=cervino_spettralo3_trig
#SBATCH --output=/scratch/jan.beutel/slurm_logs/cervino_spettralo3_trig_%A_%a.log
#SBATCH --error=/scratch/jan.beutel/slurm_logs/cervino_spettralo3_trig_%A_%a.err
###SBATCH --array=0-32%4        # 11 years × 3 channels = 33 tasks
#SBATCH --array=0-23%6        # 8 years × 3 channels = 24 tasks
# #SBATCH --time=2:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=2
# #SBATCH --gres=gpu:1
#SBATCH --partition=IFIall

#IFIgpuL40S

# Define ranges
#YEARS=(2015 2016 2017 2018 2019 2020 2021 2022 2023 2024 2025)
#STATION=MH36
YEARS=(2018 2019 2020 2021 2022 2023 2024 2025)
STATION=MH44
CHANNELS=(EHE.D EHN.D EHZ.D)

# Map task ID to year/channel
YEAR_IDX=$(( SLURM_ARRAY_TASK_ID / ${#CHANNELS[@]} ))
CHAN_IDX=$(( SLURM_ARRAY_TASK_ID % ${#CHANNELS[@]} ))

YEAR=${YEARS[$YEAR_IDX]}
CHANNEL=${CHANNELS[$CHAN_IDX]}

echo "SLURM running YEAR=$YEAR, STATION=$STATION, CHANNEL=$CHANNEL"

# Run MATLAB with year + channel arguments
matlab -batch "CERVINO_spettralo3_trig($YEAR, '$STATION', '$CHANNEL')"

echo "done."