#!/bin/bash
#SBATCH --job-name=cervino_plot_spettralo3_events
#SBATCH --output=/scratch/jan.beutel/slurm_logs/cervino_spettralo3_trig_%A_%a.log
#SBATCH --error=/scratch/jan.beutel/slurm_logs/cervino_spettralo3_trig_%A_%a.err
# #SBATCH --time=2:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=2
#SBATCH --partition=IFIall

# --- Configuration ---
CHANNELS=(EHE.D EHN.D EHZ.D)
# MH36, MH38: 2015–2025
YEARS_MH36=(2015 2016 2017 2018 2019 2020 2021 2022 2023 2024 2025)
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
    sbatch --array=0-$((TOTAL_TASKS - 1))%8 "$0"
    exit 0
fi

# --- Get job parameters ---
IDX=$SLURM_ARRAY_TASK_ID
STATION=${STATION_LIST[$IDX]}
YEAR=${YEAR_LIST[$IDX]}
CHANNEL=${CHANNEL_LIST[$IDX]}

echo "SLURM running STATION=$STATION, YEAR=$YEAR, CHANNEL=$CHANNEL"


# Run MATLAB with year + channel arguments
matlab -batch "CERVINO_plot_spettralo3_events($YEAR, '$STATION', '$CHANNEL')"

echo "done."