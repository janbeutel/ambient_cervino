function CERVINO_import_mseed(year, station, channel)
    % year comes in as a number (e.g. 2018)
    year = string(year);

    % clear all
    close all
    clearvars -except year station channel  % don’t clear year

    % year    = "2018";
    network = "1I";
    % station = "MH44";
    % channel = "EHE.D";
    
    addpath("mseed-lib");

    fprintf("Buon giorno. We are now converting hourly miniseed files...\n")

    sens = 800; % V·s/m

    % channels = ["EHE.D", "EHN.D", "EHZ.D"];
    % for j = 1:length(channels)
    %     channel = channels(j);
    
    % data_directory = "/mnt/ifi/nes/research/geophones/binaries/" + network + "/" + station + "/" + year + "/" + channel;
    data_directory = "../../../binaries/" + network + "/" + station + "/" + year + "/" + channel;

    % List all .miniseed files in the directory
    filelist = dir(fullfile(data_directory, '*.miniseed'));

    fprintf('Processing year %s for network %s station %s channel: %s\n', year, network, station, channel);
    fprintf('Number of hourly files: %d\n\n', length(filelist));
    
    savedir = "../data/" + network + "/"+ station + "/" + year + "/" + channel;
    if ~exist(savedir, 'dir')  % check if folder exists
        mkdir(savedir);         % create folder if it doesn't exist
    end
    
    % Define log file path
    logfile = "../data/" + network + "/" + station + "/" + year + "/" + channel + "_" + year + ".log";
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


    for i=1:length(filelist)     %%%%%%%%%%%%%%%%%%%%%%%%%%%
        filename=filelist(i).name
        a=ReadMSEEDFast(fullfile(filelist(i).folder, filelist(i).name));
        % a = rdmseed(fullfile(filelist(i).folder, filelist(i).name));

        fprintf('miniseed length %d',length(a))

        % Extract and concatenate all blocks
        % FS   = a(1).SampleRate;   % Hz
        % td   = vertcat(a.t);      % all timestamps
        % sigd = vertcat(a.d);      % all samples
        
        if length(a)==1    
            FS=a.sampleRate;
            td=a.matlabTimeVector;
            sigd=a.data;
        else
            FS=a(1).sampleRate;
            td=[a(1).matlabTimeVector; a(2).matlabTimeVector];
            sigd=[a(1).data; a(2).data];
        end
            
        x=td-td(1);
        sigd=sigd-mean(sigd); %demean
    
        savename=['1I' a(1).station(end) '_' a(1).channel(end) '_' num2str(a(1).dateTime.year) '_' sprintf('%02d',a(1).dateTime.month) '_' sprintf('%02d',a(1).dateTime.day) '_'];
    
        DVec=datevec(td);
        
        hh=find(DVec(:,5)==0 & DVec(:,6)==0);
        DVref=DVec(hh,:);
    
    for h=1:length(hh)
        hour=DVref(h,4);
        savenameh=[savename sprintf('%02d',hour)];
        
        if h<length(hh)
            sig=sigd(hh(h):(hh(h+1)-1)); %microV (10-6)
            t=td(hh(h):(hh(h+1)-1));
        elseif h==24
            sig=sigd(hh(h):end);   %microV (10-6)
            t=td(hh(h):end);
        end
    
        % gain_complessivo=640*10^6 % counts/(m/s)X CERVINO DAL DATA LOGGER
        
        % Correzione della sensibilità (da μV a m/s)
        sig=sig.*1E-6/sens;    % sensitivity correction --> sig in m/s %% sig=sig/gain_complessivo (m/s) questa è la conversione giusta considerando che i dati sono stored come counts e non microV

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
    
        clearvars -except filelist savedir sens i network station year channel channels fid
    end

    % Close logfile
    fclose(fid);
end
