# Ambient Cervino
Matterhorn seismic analysis


## Getting started
Tools for processing Matterhorn seismic data. Typically one year, per station, per channel. Intermediate formats in ```.mat``` files.

### Toolflow
- ```CERVINO_import_mseed_orari```: Importing waveform data into matlab, stores .mat files
- ```CERVINO_fft```: Calculate yearly FFTs
- ```CERVINO_event_detection```: Detection of events, storing of event waveforms as per-event .mat files
- ```CERVINO_plot_spettralo3_events```: Plot of event spectral properties

- ```CERVINO_classificazione```: Calculate Hilbert based classification of events
- ```CERVINO_plot_hv```: Plot stuff
- ```CERVINO_plot_hilbert```: Plot stuff

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
(c) 2025
Valeria Strallo
Chiara Colombero
Samuel Weber
Jan Beutel

## License
For open source projects, say how it is licensed.

