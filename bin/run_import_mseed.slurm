#!/bin/bash
#SBATCH --job-name=cervino_import_mseed
#SBATCH --output=/scratch/jan.beutel/slurm_logs/cervino_import_mseed_%A_%a.log
#SBATCH --error=/scratch/jan.beutel/slurm_logs/cervino_import_mseed_%A_%a.err
##SBATCH --array=0-32%6         # 11 years × 3 channels = 33 tasks
#SBATCH --array=0-5%9          # 8 years × 3 channels = 24 tasks
###SBATCH --time=12:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=2
#SBATCH --partition=IFIall

NETWORK=1I
CHANNELS=(EHE.D EHN.D EHZ.D)

# Define station-year-location mapping
declare -A STATION_YEARS
declare -A STATION_LOCATIONS

STATION_YEARS["MH36"]="2015 2016 2017 2018 2019 2020 2021 2022 2023 2024 2025"
STATION_LOCATIONS["MH36"]="A"

STATION_YEARS["MH38"]="2015 2016 2017 2018 2019"
STATION_LOCATIONS["MH38"]="A"

STATION_YEARS["MH44"]="2018 2019 2020 2021 2022 2023 2024 2025"
STATION_LOCATIONS["MH44"]="A"

STATION_YEARS["MH48"]="2019 2020 2021 2022 2023 2024 2025"
STATION_LOCATIONS["MH48"]="A B"

STATION_YEARS["MH54"]="2019 2020"
STATION_LOCATIONS["MH54"]="A"

# Flatten combinations into an array
COMBINATIONS=()
for STATION in "${!STATION_YEARS[@]}"; do
  YEARS=(${STATION_YEARS[$STATION]})
  LOCS=(${STATION_LOCATIONS[$STATION]})
  for YEAR in "${YEARS[@]}"; do
    for LOC in "${LOCS[@]}"; do
      for CHAN in "${CHANNELS[@]}"; do
        COMBINATIONS+=("$YEAR|$STATION|$LOC|$CHAN")
      done
    done
  done
done

# Get current combination
COMBO="${COMBINATIONS[$SLURM_ARRAY_TASK_ID]}"
IFS='|' read -r YEAR STATION LOCATION CHANNEL <<< "$COMBO"

echo "SLURM running import for YEAR=$YEAR, STATION=$STATION, LOCATION=$LOCATION, CHANNEL=$CHANNEL"

# Run MATLAB
matlab -batch "CERVINO_import_mseed_orari('$NETWORK', $YEAR, '$STATION', '$LOCATION', '$CHANNEL')"

echo "done."

