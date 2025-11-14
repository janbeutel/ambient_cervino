% Matlab Miniseed data importer with denoising
%
% (c) University of Innsbruck Valeria Strallo, Jan Beutel, 2025
%
%

function CERVINO_import_mseed_orari(network, year, station, location, channel)
    % year comes in as a number (e.g. 2018)
    year = string(year);

    % clear all
    close all
    clearvars -except network year station location channel  % don’t clear year

    % year    = "2018";
    % network = "1I";
    % station = "MH44";
    % channel = "EHE.D";
    % location = "A";
    
    addpath("mseed-lib");

    fprintf("Buon giorno. We are now converting hourly miniseed files...\n")

    sens = 800; % V·s/m

    % channels = ["EHE.D", "EHN.D", "EHZ.D"];
    % for j = 1:length(channels)
    %     channel = channels(j);
    
    % data_directory = "/mnt/ifi/nes/research/geophones/binaries/" + network + "/" + station + "/" + year + "/" + channel;
    data_directory = "../../../binaries/" + network + "/" + station + "/" + year + "/" + channel;

    % List all .miniseed files in the directory
    filelist = dir(fullfile(data_directory, network + "." + station + "." + location + "." + channel + '*.miniseed'));

    fprintf('Processing year %s for network %s station %s location %s channel: %s\n', year, network, station, location, channel);
    fprintf('Number of hourly files: %d\n\n', length(filelist));
    
    savedir = "../data/" + network + "/"+ station + "/" + year + "/" + channel;
    if ~exist(savedir, 'dir')  % check if folder exists
        mkdir(savedir);         % create folder if it doesn't exist
    end
    
    % Define log file path
    logfile = "../data/" + network + "/" + station + "/" + year + "/" + location + "." + channel + "_" + year + ".log";
    % Remove logfile if it exists
    if exist(logfile, "file")
        delete(logfile);
    end
    % Open in append mode ("a")
    fid = fopen(logfile, "a");
    if fid == -1
        error("Could not open log file.");
    end
    % Write header logfile
    fprintf(fid, "miniseedfilename,numsamples,samplerate,matfilename\n");

    for i = 1:length(filelist)
        filename = filelist(i).name;
    
        a = rdmseed(fullfile(filelist(i).folder, filelist(i).name));
        % disp(a(1))
    
        % Extract and concatenate all blocks
        FS   = a(1).SampleRate;   % Hz
        td   = vertcat(a.t);      % all timestamps
        sigd = vertcat(a.d);      % all samples
    
        % Rimozione media
        sigd = sigd - mean(sigd);
        
        % Tempo relativo
        % t = td - td(1);
        t = td;
    
        % Correzione della sensibilità (da μV a m/s)
        sig = sigd * 1e-6 / sens;
    
        fprintf("%s numsamples %d sample rate %d\n", filename, length(sig), FS)
    
        % Verifica se il file contiene 1 ora di dati
        if length(sig) == FS * 3600;
            savename = erase(filename, '.miniseed') + ".mat";
            save(savedir + "/" + savename, 'sig', 't', 'FS');

            % Write line to logfile
            fprintf(fid, "%s,%d,%d,%s\n", filename, length(sig), FS, savename);
        else
            warning('File %s non contiene esattamente 1 ora di dati %d != %d samples. Skippato.', filename, length(sig), FS * 3600);
            % Write line to logfile
            fprintf(fid, "%s,%d,%d,\n", filename, length(sig), FS);
        end
    
        clearvars -except filelist savedir sens i network station year location channel channels fid
    end

    % Close logfile
    fclose(fid);
end