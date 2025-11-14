function CERVINO_fft(network, year, station, location, channel)
    addpath(genpath('lib'));

    % year comes in as a number (e.g. 2018)
    year = string(year);

    close all
    clearvars -except network year station location channel  % don’t clear year

    % year    = "2018";
    % network = "1I";
    % station = "MH44";
    % channel = "EHE.D";

    data_directory = "../data/" + network + "/" + station + "/" + year + "/" + channel;
    % List all .miniseed files in the directory
    filelist = dir(fullfile(data_directory, network + "." + station + "." + location + "." + channel + '*.mat'));
    TOT = size(filelist,1);

    savedir = "../results/fft/" + network + "/"+ station + "/" + year;
        if ~exist(savedir, 'dir')  % check if folder exists
            mkdir(savedir);        % create folder if it doesn't exist
        end

    fprintf("Running CERVINO_fft for year %s, station %s, channel %s, %d\n", year, station, channel, TOT)

    FFT_all=zeros(500,TOT);

    for ite= 1:1:TOT
        fprintf("File %d, station %s, channel %s, %s\n", ite, station, channel, filelist(ite).name)

        % get the file name:
        folder = filelist(ite).folder;
        filename = filelist(ite).name;
        load(folder + "/" + filename); 

        DVec(ite)=t(1);

        L=length(sig);
        sig=sig-mean(sig);   %demean the signal before further computations
        
        %FS=250;
        
        % Parameters
        Fnyquist=FS/2;
        std_range=4;
        l_win=300;
        nfft=2^nextpow2(l_win*FS);
        apodisation=0.1;
        cosine_taper=tukeywin(nfft,apodisation);
        
        % Signal    
        duV=sig;
        tuV=t-t(1);
        
        % f vector
        vect_fz=(FS/2)*linspace(0,1,nfft/2+1);
        freq_band=vect_fz<=Fnyquist;
        vect_fz_band=vect_fz(freq_band);

        nb_win=floor((length(duV)-1)/(l_win*FS));     %l_win [s]
        mFFT=zeros(500,nb_win-1);                     %Initiation de la matrice mFFT
        
        for w=1:nb_win-1
            %Decoupage du signal en fenetre de l_win secondes
            signal_fenetre=duV((w-1)*l_win*FS+w:w*l_win*FS+w);
            
            %Calcul fft
            [vect_f,fftV] = func_spectre_puissance_c(signal_fenetre,FS,nfft,apodisation);
            clear signal_fenetre
            
            %Interpolation logarithmique du spectre entre Fmin et Fmax Hz sur 2000 pts (diminue le temps de calcul du lissage)
            Fmin=0.1;
            Fmax=100;
            vect_f_int=logspace(log10(Fmin),log10(Fmax),1000)';
            
            FFTV=interp1(vect_f,fftV,vect_f_int,'spline');
            clear vect_f fftV
            
            %Lissage des spectres par la technique de Konno et Ohmachi (1998) avec un parametre de lissage b=40
            bandwidth=70;
            FFTV_liss=konnoOhmachiSmoothing(FFTV,vect_f_int,bandwidth,'normalize');
            
            %Interpolation logarithmique du spectre entre Fmin et Fmax Hz avec un pas de 500 pts
            VECT_F=logspace(log10(Fmin),log10(Fmax),500)';
            FFT=interp1(vect_f_int,FFTV_liss,VECT_F,'spline');
            mFFT(:,w)=FFT;
            clear fftV
        end  
        
        cFFT=mean(mFFT,2);   
        
        FFT_all(:,ite)=cFFT;
        
        clearvars -except FFT_all VECT_F DVec filelist ite TOT VECT_F savename savedir network station year location channel

    end

    % save(savename,'FFT_all','DVec','VECT_F')

    %save CLASSE_MH44_provan 
    % save CLASS_MH44_N Hilbert Results DVec1

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    savename=["fft_" + station + "." + location + "." + channel + "." + year + ".mat"];
save(savedir + "/" + savename, "FFT_all","DVec","VECT_F")
end