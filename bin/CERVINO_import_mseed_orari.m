% Matlab Miniseed data importer with denoising
%
% (c) University of Innsbruck Valeria Strallo, Jan Beutel, 2025
%
%

fprintf("Buon giorno. We are now converting hourly miniseed files...\n")

clear all
close all

addpath("mseed-lib");
% savepath;

sens = 800; % V·s/m

network = "1I";
station = "MH44";
channel = "EHN.D";
year    = "2024";

channels = ["EHE.D", "EHN.D", "EHZ.D"];

for j = 1:length(channels)
    channel = channels(j);
    
    data_directory = "../../../binaries/" + network + "/" + station + "/" + year + "/" + channel;

    % List all .miniseed files in the directory
    filelist = dir(fullfile(data_directory, '*.miniseed'));

    fprintf('Processing year %s for network %s station %s channel: %s\n', year, network, station, channel);
    fprintf('Number of hourly files: %d\n\n', length(filelist));
    
    for i = 1:length(filelist)
        filename = filelist(i).name;
        disp(filename)   % show current file
    
        a = rdmseed(fullfile(filelist(i).folder, filelist(i).name));
        % disp(a(1))
    
        % Extract and concatenate all blocks
        FS   = a(1).SampleRate;   % Hz
        td   = vertcat(a.t);      % all timestamps
        sigd = vertcat(a.d);      % all samples
    
        % Rimozione media
        sigd = sigd - mean(sigd);
        
        % Tempo relativo
        t = td - td(1);
    
        % Correzione della sensibilità (da μV a m/s)
        sig = sigd * 1e-6 / sens;
    
        % lengthsignal = length(sig)
        % expected = FS * 3600
    
        % Verifica se il file contiene 1 ora di dati
        if length(sig) == FS * 3600;
            savename = erase(filename, '.miniseed') + ".mat";
        
            savedir = "../data/" + network + "/"+ station + "/" + year + "/" + channel;
            if ~exist(savedir, 'dir')  % check if folder exists
                mkdir(savedir);         % create folder if it doesn't exist
            end
            save(savedir + "/" + savename, 'sig', 't', 'FS');
        else
            warning('File %s non contiene esattamente 1 ora di dati %d != %d samples. Skippato.', filename, length(sig), FS * 3600);
        end
    
        clearvars -except filelist sens i network station year channel channels
    end
end