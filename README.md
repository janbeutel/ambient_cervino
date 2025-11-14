# Ambient Cervino
Matterhorn seismic analysis toolbox


## Getting started
Tools for processing Matterhorn seismic data. Processing is typically done per one year, per station, per location, per channel. Intermediate formats are stored in ```.mat``` files.

### Toolflow
- ```CERVINO_import_mseed_orari```: Importing waveform data into matlab, stores .mat files
- ```CERVINO_fft```: Calculate annual FFTs
- ```CERVINO_event_detection```: Detection of events, storing of event waveforms as per-event .mat files
- ```CERVINO_plot_spettralo3_events```: Plot of per event spectral properties

- ```CERVINO_classificazione```: Calculate Hilbert based classification of events
- ```CERVINO_plot_hv```: Plot stuff
- ```CERVINO_plot_hilbert```: Plot stuff

### Automation using SLURM
SLURM is the cluster batch processing system that is used to automate batch processing accross array jobs. Each processing step contains a *.slurm file like ```sbatch run_event_detection.slurm``` to start a SLURM array job. The SLURM queue details need to be adapted for each specific cluster instance. 

## Repository structure
- ```bin```: analysis toolbox code
- ```data```: space for input data and converted products
- ```results```: analysis products and plots
- ```tools```: further analysis tools


## ToDo
- Do we want to process partial hourly files?

## Authors and acknowledgment
(c) 2025
* Valeria Strallo - Politecnico di Torino
* Chiara Colombero - Politecnico di Torino
* Samuel Weber - SLF Davos
* Jan Beutel - University of Innsbruck


