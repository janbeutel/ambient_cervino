# Ambient Cervino
Matterhorn seismic analysis


## Getting started
Tools for processing Matterhorn seismic data. Typically one year, per station, per channel. Intermediate formats in ```.mat``` files.

### Toolflow
- ```CERVINO_import_mseed_orari```: Importing waveform data into matlab, stores .mat files
- ```CERVINO_event_detection```: Detection of events, storing of event waveforms as .mat files
- ```CERVINO_spettralo3_trig```: Plot of spectral properties

### Automation using SLURM
SLURM is the cluster batch processing system that is used to automatebatch processing accross array jobs. Typically per STATION, CHANNEL, YEAR

```sbatch run_event_detection.slurm```

## Repository structure
- ```bin```: code
- ```data```: space for input data and converted products
- ```results```: analysis products


## ToDo
- Opensource? License?
- Do we want to process partial hourly files?

## Authors and acknowledgment
Valeria Strallo
Chiara Colombero
Samuel Weber
Jan Beutel

## License
For open source projects, say how it is licensed.

