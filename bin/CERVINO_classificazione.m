function CERVINO_classificazione(year, station, channel)
   % year comes in as a number (e.g. 2018)
   year = string(year);

   close all
   clearvars -except year station channel  % don’t clear year

   % year    = "2018";
   network = "1I";
   % station = "MH44";
   % channel = "EHE.D";

    data_directory = "../results/events/" + network + "/" + station + "/" + year + "/" + channel;
    % List all .miniseed files in the directory
    FileList = dir(fullfile(data_directory, '*.mat'));
    TOT = size(FileList,1);
    
    savedir = "../results/classification/" + network + "/"+ station + "/" + year;
    if ~exist(savedir, 'dir')  % check if folder exists
        mkdir(savedir);         % create folder if it doesn't exist
    end

    fprintf("Running CERVINO_classificazione for year %s, station %s, channel %s, %d\n", year, station, channel, TOT)

    for ite= 1:1:TOT
        fprintf("Event %d, station %s, channel %s, %s\n", ite, station, channel, FileList(ite).name)

        % get the file name:
        folder = FileList(ite).folder;
        filename = FileList(ite).name;
        load(folder + "/" + filename) 

        L=length(SIG);
        SR=1/FS;

        year=str2num(filename(17:20));
        month=str2num(filename(21:22));
        day=str2num(filename(23:24));
        hour=str2num(filename(26:27));
        mm=floor(sigseconds/60);
        ss=floor(sigseconds-(mm*60));
            
        DVec1(ite,:)=[year month day hour mm ss];
        % DVec1(ite,:)=[year month day hour];

        t=(X-X(1))*24*3600;
        sig=SIG-mean(SIG);

        %% MAX AMPLITUDE OF THE SIGNAL (hilbert 10)
        A_N=max(SIG(:,1));
        Hilbert(ite,10)=A_N;

        %%% 1. Ratio of the Maximum Amplitude to the Mean of the Envelope P1 -
        %%% hilbert 1
        y=hilbert(sig);
        envelope=abs(y);

        Amedia=mean(envelope);
        P1=max(sig)/Amedia;

        %%% 2. Kurtosis of the envelope P2
        %%% hilbert 2
        P2=kurtosis(sig);

        Hilbert(ite,1)=P1;
        Hilbert(ite,2)=P2;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Duration

        datamax2=max(envelope);
        timax=find(envelope==datamax2);

        %Set the threshold (15% A max)
        threshold=datamax2*0.15;

        % Trova gli indici in cui il segnale supera la soglia
        above_threshold = find(envelope >= threshold);

        % Calcola la Bracketed Duration (da primo a ultimo campione sopra soglia) P3
        if isempty(above_threshold)
            bracketed_dur = 0;
            uniform_dur = 0;
        else
            bracketed_dur = (above_threshold(end) - above_threshold(1)) / FS;

            % Calcola la Uniform Duration (somma durata intervalli sopra soglia) P4
            uniform_dur = numel(above_threshold) / FS;
        end

        P3=bracketed_dur;
        P3bis=uniform_dur;

        % Lenght of the event file
        P3ter=t(end);

        % Increasing/Decreasing time ratio
        a=t(above_threshold(1));
        b=t(timax);
        c=t(above_threshold(end));

        t_increasing=b-a;
        t_decreasing=c-b;
        IncDec=t_increasing/t_decreasing;

        P4=IncDec;

        Hilbert(ite,3)=P3;
        Hilbert(ite,4)=P3bis;
        Hilbert(ite,5)=P3ter;
        Hilbert(ite,6)=P4;


        % Amplitude/duration ratio hilbert 11
        Hilbert(ite,11)=A_N/P3bis;

        %%% 5. Energy in the 10-30 Hz frequency Band (1-->5 hilbert et al. 2014
        %%% Automated Indentification, location, and volume estimation of rockfalls
        %%% at Piton de la Fournaise volcano)

        t=0 : 1/FS : 1-(1/FS);
        freq = FS*(0:(L/2))/L;

        Results(ite,2)=sigseconds;

        SR=1/FS;
        freq = 0 : (1/(length(sig)*SR)) : 0.5/(SR);

        sp1 = abs(fft(sig));
        sp1 = sp1(1:length(freq));
        [M1, Indices1] = max(sp1(2:end));
        fpicco_1=freq(Indices1);

        P6= fpicco_1;
        Hilbert(ite,9)=P6;

        % Width of the spectrum above 50% of the peak frequency P5bis hilbert 8

        %Set the threshold (50% A max)
        threshold2=M1*0.80;

        % Trova gli indici dove lo spettro supera la soglia
        above_threshold = sp1 >= threshold2;

        % Identifica gli intervalli di frequenza sopra la soglia
        df = diff(freq); % Passo di frequenza tra campioni
        bandwidths = [];

        % Identifica gli intervalli di frequenza sopra la soglia
        start_idx = find(diff([0, above_threshold']) == 1); % Indici di inizio
        end_idx = find(diff([above_threshold', 0]) == -1); % Indici di fine

        % Calcola la lunghezza di banda per ogni intervallo
        for b = 1:length(start_idx)
            bandwidths(b) = freq(end_idx(b)) - freq(start_idx(b));
        end

        % Somma le lunghezze di banda
        total_bandwidth = sum(bandwidths);

        P5bis=total_bandwidth;
        Hilbert(ite,8)=P5bis;

        clearvars -except k FileList TOT savedir Ampiezze Ampiezze_ST Results Amax_ST Hilbert Parametri DVec1 FS names network station year channel
        %clearvars -except DVec1 k FileList1 TOT 
    end 


    %save CLASSE_MH44_provan 
    % save CLASS_MH44_N Hilbert Results DVec1

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    savename=["classification_" + station + "_" + channel + "_" + year + ".mat"];
    save(savedir + "/" + savename, "Hilbert", "Results", "DVec1")
end