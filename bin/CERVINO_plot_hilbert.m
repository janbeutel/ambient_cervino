function CERVINO_plot_hilbert(network, year, station, location)
    addpath(genpath('lib'));

    % year comes in as a number (e.g. 2018)
    year = string(year);

    close all
    clearvars -except network year station location channel  % don’t clear year

    % year    = "2018";
    % network = "1I";
    % station = "MH44";
    % channel = "EHE.D";

    data_directory = "../results/fft/" + network + "/" + station + "/" + year + "/";
    % List all .miniseed files in the directory
    filelist = dir(fullfile(data_directory, "fft_" + station + "." + location + '*.mat'));
    TOT = size(filelist,1);

    savedir = "../results/plots/";
    if ~exist(savedir, 'dir')  % check if folder exists
        mkdir(savedir);         % create folder if it doesn't exist
    end

    fprintf("Running CERVINO_plot_hilbert for network %s, year %s, station %s, location %s, num files %d\n", network, year, station, location, TOT)


    channel = "EHE.D";
    filename=["fft_" + station + "." + location + "." + channel + "." + year + ".mat"];
    load(data_directory + filename) 
    FFT_E=FFT_all;
    DV_E=DVec;

    channel = "EHN.D";
    filename=["fft_" + station + "." + location + "." + channel + "." + year + ".mat"];
    load(data_directory + filename) 
    FFT_N=FFT_all;
    DV_N=DVec;

    channel = "EHZ.D";
    filename=["fft_" + station + "." + location + "." + channel + "." + year + ".mat"];
    load(data_directory + filename) 
    FFT_Z=FFT_all;
    DV_Z=DVec;


    DVc=DV_E(1):1/24:DV_E(end);   % complete vector
    FFT_E_c=zeros(length(VECT_F),length(DVc));
    FFT_N_c=zeros(length(VECT_F),length(DVc));
    FFT_Z_c=zeros(length(VECT_F),length(DVc));

        for i=1:length(DV_E)
            e=find((abs(DVc-DV_E(i)))==min(abs(DVc-DV_E(i))));
            FFT_E_c(:,e)=FFT_E(:,i);
        end

        for i=1:length(DV_N)
            e=find((abs(DVc-DV_N(i)))==min(abs(DVc-DV_N(i))));
            FFT_N_c(:,e)=FFT_N(:,i);
        end
        
        for i=1:length(DV_Z)
            e=find((abs(DVc-DV_Z(i)))==min(abs(DVc-DV_Z(i))));
            FFT_Z_c(:,e)=FFT_Z(:,i);
        end
        
    %Parameters
    type_HV='THE';  %Three differents combinations : 'SA'=sqared average  or  'THE'=total horizontal energy  or  'DE'=directional energy or 'E'=Chamousset2_before_2014/05/27
    theta=0;  %degree

    %Limites axes
    Pmin=0.1;   %HV min axe
    Pmax=10;      %HV max axe
    Res=100;        %Resolution of PDF
    THETA=theta*pi/180;
    Theta=num2str(theta);

    fmin=0.1; %Hz
    fmax=100; %Hz

    %Bins definition
    %vecbinf=logspace(-3,2,100); %vecteur espace regulierement en log, entre 10^-3 et 100 Hz, avec 10 points par decade
    % vecbinf=(0.5:0.1:100); %vecteur espace regulierement en linéaire de 0.1 Hz, entre 0.01 et 100 Hz
    vecbinf=logspace(log10(fmin),log10(fmax),200);

    vecbinf=vecbinf(:);
    vecbinpow=logspace(log10(0.0625),log10(512),Res); %vecteur espace regulierement en lineaire, entre -300 et 0 dB
    vecbinpow=vecbinpow(:);

    %Load vecteur f
    %     vect_f=logspace(log10(fmin),log10(fmax),500)';
    vect_f=VECT_F;

    %Initialisation matrice PDF
    Super_PDF_mat=zeros(length(vecbinpow),length(vecbinf));
    Stats_HV_mat=zeros(size(FFT_E_c));
            
                    %Selection de la composante verticale, nord et est
                    Super_FFT_Z=FFT_Z_c;
                    Super_fft_N=FFT_N_c;
                    Super_fft_E=FFT_E_c;
                    
                        for ii=1:length(Super_FFT_Z(1,:));
                            SV=Super_FFT_Z(:,ii);
                            SN=Super_fft_N(:,ii);
                            SE=Super_fft_E(:,ii);
                            
                            % Spectres horizontal
                            switch type_HV
                                case 'SA';   SH=sqrt((SN.^2+SE.^2)/2);
                                case 'THE';  SH=sqrt(SN.^2+SE.^2);
                                case 'DE';
                                    if isnan(theta)
                                        disp('Warning theta not defined');
                                    else
                                        SH=(SN.*cos(THETA))+(SE.*sin(THETA));
                                    end
                                case 'E'
                                    SH=SE;
                            end
                            
                            if sum(SV)~=0 && sum(SH)~=0
                                HV=SH./SV;
                                HV=HV(:);                 
                                if max(HV)<=512 && min(HV)>=0.0625
    %                                 Stats_HV_mat=[Stats_HV_mat HV];
                                [Matrix_Counts] = func_hist2d_noc(vect_f,HV,vecbinf,vecbinpow);
                                Super_PDF_mat=Super_PDF_mat+Matrix_Counts;
                                Stats_HV_mat(:,ii)=HV;
                                end
                            end
                        end

    %Normalisation probas
    N_occurences_tot=sum(Super_PDF_mat,1);

    %Normalisation
    for ind_f=1:length(vecbinf)
        Super_PDF_mat(:,ind_f)=Super_PDF_mat(:,ind_f)/N_occurences_tot(ind_f);
    end

    %Stats
    a=find(Stats_HV_mat(100,:)==0);
    mat=Stats_HV_mat;                             
    mat(:,a)=[];
    DV_short=DVc;
    DV_short(a)=[];


    % channel = "EHN.D";
    channels = ["EHE.D", "EHN.D", "EHZ.D"];

    for channel = channels

        % Example: your existing code adapted to this loop
        disp("Plotting channel: " + channel);

        data_directory = "../results/classification/" + network + "/" + station + "/" + year + "/";
        % List all .miniseed files in the directory
        filelist = dir(fullfile(data_directory, network + "." + station + "." + location + "." + channel + '*.mat'));
        TOT = size(filelist,1);

        filename=["classification_" + station + "." + location + "." + channel + "." + year + ".mat"];
        load(data_directory + filename) 

        load('../metadata/fieldwork_dates.mat');

        start_fw=datenum(start_fw);
        end_fw=datenum(end_fw);

        data=Hilbert;
        Dv=datenum(DVec1);

        minDv=min(Dv);
        maxDv=max(Dv);
        
        % minDv_raw = min(Dv);
        % maxDv_raw = max(Dv);

        % % Extract year from minDv
        % minyear = datevec(minDv_raw);

        % % Reset minDv and maxDv to full calendar year
        % minPlot = datenum(minyear(1), 1, 1);      % January 1st
        % maxPlot = datenum(minyear(1), 12, 31);    % December 31st

        data_f=data;
        Dv_f=Dv;
        a=find(Hilbert(:,9)>35);
        data_f(a,:)=[];
        Dv_f(a)=[];

        figure;
        tiledlayout(3,1, 'TileSpacing', 'compact', 'Padding', 'compact');
        % Increase figure height for better visibility
        set(gcf, 'Position', [100 100 900 900])  % [left bottom width height]

        % ---- Plot 1 ----
        nexttile
        plot(Dv, data(:,9),'.')
        datetick('x','yyyy-mm-dd', 'keeplimits')
        xlim([minDv maxDv])
        % xlabel('Date (mm/yy)')
        ylabel('F(Hz)')
        hold on
        for i = 1:length(start_fw)
            idx = Dv >= start_fw(i) & Dv <= end_fw(i);
            plot(Dv(idx), data(idx,9), 'r.', 'MarkerSize', 12)
        end

        % ---- Plot 2 ----
        nexttile
        plot(Dv_f, data_f(:,9),'.')
        datetick('x','mm/yy', 'keeplimits')
        xlim([minDv maxDv])
        % xlabel('Date (mm/yy)')
        ylabel('F(Hz)')
        hold on
        % hold off

        % ---- Plot 3 ----
        nexttile
        imagesc(DVc,VECT_F,Stats_HV_mat)
        set(gca,'Yscale','log')
        set(gca,'Ydir','normal')
        xlim([minDv maxDv])
        ylim([5 50])
        datetick('x','mm/yy','keeplimits')
        % xlabel('Date (mm/yy)')
        ylabel('Frequency (Hz)')
        load('lib/FFT_cmap.mat')
        colormap(cmap);
        c=colorbar;
        c.Label.String='HV';
        caxis([0 10])
        hold on
        set(gca,'layer','top')
        % box on
        % sdf('20_15')

        savename = "hilbert1_" + station + "." + location + "." + channel + "." + year;
        saveas(gcf, fullfile(savedir, savename + ".jpg"))
        %  saveas(gcf,[savename '_norm.fig'])  


        figure;
        tiledlayout(3,1, 'TileSpacing', 'compact', 'Padding', 'compact');
        % Increase figure height for better visibility
        set(gcf, 'Position', [100 100 900 900])  % [left bottom width height]

        % ---- Plot 1 ----
        nexttile
        plot(Dv, data(:,4),'.')
        datetick('x','mm/yy', 'keeplimits')
        xlim([minDv maxDv])
        % xlabel('Date (mm/yy)')
        ylabel('Uniform duration (s)')
        hold on
        for i = 1:length(start_fw)
            idx = Dv >= start_fw(i) & Dv <= end_fw(i);
            plot(Dv(idx), data(idx,4), 'r.', 'MarkerSize', 12)
        end

        % ---- Plot 2 ----
        nexttile
        plot(Dv, data(:,10),'.')
        datetick('x','mm/yy', 'keeplimits')
        xlim([minDv maxDv])
        ylim([0,0.0002])
        % xlabel('Date (mm/yy)')
        ylabel('Max amplitude (dB)')
        hold on
        for i = 1:length(start_fw)
            idx = Dv >= start_fw(i) & Dv <= end_fw(i);
            plot(Dv(idx), data(idx,10), 'r.', 'MarkerSize', 12)
        end

        % ---- Plot 3 ----
        nexttile
        plot(Dv, data(:,1),'.')
        datetick('x','mm/yy', 'keeplimits')
        xlim([minDv maxDv])
        % xlabel('Date (mm/yy)')
        ylabel('Amax/Amean')
        hold on
        for i = 1:length(start_fw)
            idx = Dv >= start_fw(i) & Dv <= end_fw(i);
            plot(Dv(idx), data(idx,1), 'r.', 'MarkerSize', 12)
        end
        % sdf('20_15')

        savename = "hilbert2_" + station + "." + location + "." + channel + "." + year;
        saveas(gcf, fullfile(savedir, savename + ".jpg"))
        %  saveas(gcf,[savename '_norm.fig'])  
    end
end


