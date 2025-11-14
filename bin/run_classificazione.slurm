#!/bin/bash
#SBATCH --job-name=cervino_classificazione
#SBATCH --output=/scratch/jan.beutel/slurm_logs/cervino_classificazione_%A_%a.log
#SBATCH --error=/scratch/jan.beutel/slurm_logs/cervino_classificazione_%A_%a.err
#SBATCH --time=2:00:00
#SBATCH --mem=25G
#SBATCH --cpus-per-task=4
#SBATCH --partition=IFIgpuL40S

# --- Configuration ---
NETWORK=1I
CHANNELS=(EHE.D EHN.D EHZ.D)
YEARS_MH36=(2015 2016 2017 2018 2019 2020 2021 2022 2023 2024 2025)
YEARS_MH38=(2015 2016 2017 2018 2019)
YEARS_MH44=(2018 2019 2020 2021 2022 2023 2024 2025)
YEARS_MH48=(2019 2020 2021 2022 2023 2024 2025)
YEARS_MH54=(2019 2020)

# Define locations for each station
declare -A STATION_LOCATIONS
STATION_LOCATIONS[MH36]="A"
STATION_LOCATIONS[MH38]="A"
STATION_LOCATIONS[MH44]="A"
STATION_LOCATIONS[MH48]="A B"
STATION_LOCATIONS[MH54]="A"

# --- Build job mapping ---
declare -a STATION_LIST
declare -a YEAR_LIST
declare -a CHANNEL_LIST
declare -a LOCATION_LIST

for STATION in MH36 MH38 MH44 MH48 MH54; do
    eval YEARS_ARRAY=( "\${YEARS_${STATION}[@]}" )
    LOCATIONS=( ${STATION_LOCATIONS[$STATION]} )  # split into array
    for YEAR in "${YEARS_ARRAY[@]}"; do
        for CHANNEL in "${CHANNELS[@]}"; do
            for LOCATION in "${LOCATIONS[@]}"; do
                STATION_LIST+=("$STATION")
                YEAR_LIST+=("$YEAR")
                CHANNEL_LIST+=("$CHANNEL")
                LOCATION_LIST+=("$LOCATION")
            done
        done
    done
done

TOTAL_TASKS=${#STATION_LIST[@]}

# --- SLURM dynamic setup ---
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
LOCATION=${LOCATION_LIST[$IDX]}

echo "SLURM running classification for YEAR=$YEAR, STATION=$STATION, LOCATION=$LOCATION, CHANNEL=$CHANNEL"

# Run MATLAB with year + channel arguments
matlab -batch "CERVINO_classificazione('$NETWORK', $YEAR, '$STATION', '$LOCATION', '$CHANNEL')"

echo "done."