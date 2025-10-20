#!/bin/bash
#SBATCH --job-name=cervino_event_detection
#SBATCH --output=/scratch/jan.beutel/slurm_logs/cervino_event_detection_%A_%a.log
#SBATCH --error=/scratch/jan.beutel/slurm_logs/cervino_event_detection_%A_%a.err
# SBATCH --time=12:00:00
#SBATCH --mem=25G
#SBATCH --cpus-per-task=4
#SBATCH --partition=IFIall
# SBATCH --array=0-??%8   # will be set dynamically below

# --- Configuration ---
CHANNELS=(EHE.D EHN.D EHZ.D)
# MH36, MH38: 2015–2025
YEARS_MH36=(2015 2016 2017 2018 2019 2020 2021 2022 2023 2024 2025)
#YEARS_MH36=(2021 2022 2023 2024 2025)
YEARS_MH38=(2015 2016 2017 2018 2019)
# MH44: only from 2018
YEARS_MH44=(2018 2019 2020 2021 2022 2023 2024 2025)

# --- Build job mapping ---
declare -a STATION_LIST
declare -a YEAR_LIST
declare -a CHANNEL_LIST

for STATION in MH44 MH36 MH38; do
    # pick station-specific years
    eval YEARS_ARRAY=( "\${YEARS_${STATION}[@]}" )
    for YEAR in "${YEARS_ARRAY[@]}"; do
        for CHANNEL in "${CHANNELS[@]}"; do
            STATION_LIST+=("$STATION")
            YEAR_LIST+=("$YEAR")
            CHANNEL_LIST+=("$CHANNEL")
        done
    done
done

TOTAL_TASKS=${#STATION_LIST[@]}

# --- SLURM dynamic setup ---
# SLURM does not parse variables in #SBATCH lines, so we re-submit dynamically if not yet running under SLURM
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
    echo "Submitting ${TOTAL_TASKS} jobs (max 8 concurrent)..."
    sbatch --array=0-$((TOTAL_TASKS - 1))%10 "$0"
    exit 0
fi

# --- Get job parameters ---
IDX=$SLURM_ARRAY_TASK_ID
STATION=${STATION_LIST[$IDX]}
YEAR=${YEAR_LIST[$IDX]}
CHANNEL=${CHANNEL_LIST[$IDX]}

echo "SLURM running STATION=$STATION, YEAR=$YEAR, CHANNEL=$CHANNEL"

# --- Run MATLAB task ---
matlab -batch "CERVINO_event_detection($YEAR, '$STATION', '$CHANNEL')"

echo "done."
