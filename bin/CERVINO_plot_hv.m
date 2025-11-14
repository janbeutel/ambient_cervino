function CERVINO_plot_hv(network, year, station, location)
    addpath(genpath('lib'));

    % year comes in as a number (e.g. 2018)
    year = string(year);

    close all
    clearvars -except network year station location channel  % don’t clear year

    % year    = "2018";
    network = "1I";
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

    fprintf("Running CERVINO_hv for network %s, year %s, station %s, location %s, num files %d\n", network, year, station, location, TOT)


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
    FFT_V=FFT_all;
    DV_V=DVec;


    DVc=DV_E(1):1/24:DV_E(end);   % complete vector
    FFT_E_c=zeros(length(VECT_F),length(DVc));
    FFT_N_c=zeros(length(VECT_F),length(DVc));
    FFT_V_c=zeros(length(VECT_F),length(DVc));

    for i=1:length(DV_E)
        e=find((abs(DVc-DV_E(i)))==min(abs(DVc-DV_E(i))));
        FFT_E_c(:,e)=FFT_E(:,i);
    end

    for i=1:length(DV_N)
        e=find((abs(DVc-DV_N(i)))==min(abs(DVc-DV_N(i))));
        FFT_N_c(:,e)=FFT_N(:,i);
    end

    for i=1:length(DV_V)
        e=find((abs(DVc-DV_V(i)))==min(abs(DVc-DV_V(i))));
        FFT_V_c(:,e)=FFT_V(:,i);
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
    % vecbinf=(0.5:0.1:100); %vecteur espace regulierement en lin�aire de 0.1 Hz, entre 0.01 et 100 Hz
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
    Super_fft_V=FFT_V_c;
    Super_fft_N=FFT_N_c;
    Super_fft_E=FFT_E_c;
                        
    for ii=1:length(Super_fft_V(1,:));
        SV=Super_fft_V(:,ii);
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
        
    %   Stats_HV_mat=Stats_HV_mat(:,2:end); %vire la colonne initialisation a zero
    %Mean
    mean_pow = mean(mat,2);
    mean_pow=interp1(vect_f,mean_pow,vecbinf,'linear');
    %Min
    %     min_pow = min(Stats_FFT_mat,[],2);
    %     min_pow=interp1(vect_f,min_pow,vecbinf,'linear');
    %     %Max
    %     max_pow = max(Stats_FFT_mat,[],2);
    %     max_pow=interp1(vect_f,max_pow,vecbinf,'linear');
    %Median
    %     med_pow = median(Stats_HV_mat,2);
    %     med_pow=interp1(vect_f,med_pow,vecbinf,'linear');
    %90th percentile
    nineth_pow = prctile(mat,90,2);
    nineth_pow=interp1(vect_f,nineth_pow,vecbinf,'linear');
    %10th percentile
    th_pow = prctile(mat,10,2);
    th_pow=interp1(vect_f,th_pow,vecbinf,'linear');
    %Standard deviation
    std_pow=std(mat,[],2);
    std_pow=interp1(vect_f,std_pow,vecbinf,'linear');
    min_std = mean_pow - std_pow;
    max_std = mean_pow + std_pow;
        
    %     %save
    %     clear HV
    %     HV(:,1)=mean_pow;
    %     HV(:,2)=nineth_pow;
    %     HV(:,3)=th_pow;
    %     name_hv=['hv_' Site '_' Capteur_traite '_type_HV_' type_HV '_theta_' Theta '_du_' Jour_debut '_au_' Jour_fin '_fmin_' num2str(fmin) '_fmax_' num2str(fmax) '_Hz'];
    %     nom_fichier=([P '/data_periodogram/h_hh_hv/data/' name_hv '.txt']);
    %     save(nom_fichier,'HV','-ascii','-double');
    %     
    %     Ph='../../Chamousset2/data_periodogram/h_hh_hv/data/';
    %     HV1_b=load([Ph 'hv_Chamousset2_2_type_HV_DE_theta_55_du_20120701_au_20120701_fmin_0.5_fmax_100_Hz.txt']);
    %     HV1_c=load([Ph 'hv_Chamousset2_2_type_HV_DE_theta_55_du_20120901_au_20120901_fmin_0.5_fmax_100_Hz.txt']);
    %     
        
    % Figure
    figure(1)
    hold on
    H0=pcolor(vecbinf,vecbinpow,Super_PDF_mat); hold on;
    shading flat;
    view([0,90]);
    caxis([0 max(max(Super_PDF_mat))]);
    CM=flipud(colormap(hot));
    Array_color=(5:size(CM,1));
    colormap([CM(1,:);CM(Array_color,:)]);
    caxis([0 0.5]);
    cc = colorbar;
    ylabel(cc, 'Probability (-)');
    %     text(0.7,64,[Jour_debut ' to ' Jour_fin],'HorizontalAlignment','left','BackgroundColor',[0.9 0.9 0.9],'FontSize',14,'color','red','FontWeight','bold')

    %Curves
    H3 = plot(vecbinf,mean_pow,'-k','Linewidth',2);zoom on; grid on; hold on;  
    %     plot(vecbinf,HV1_b(:,1),'b','Linewidth',2);zoom on; grid on; hold on;
    %     plot(vecbinf,HV1_c(:,1),'m','Linewidth',2);zoom on; grid on; hold on;
    %     H4 = plot(vecbinf,min_pow,'-r','Linewidth',1);zoom on; grid on; hold on;
    %     H5 = plot(vecbinf,max_pow,'-k','Linewidth',1);zoom on; grid on; hold on;
    %     H6 = plot(vecbinf,med_pow,'-b','Linewidth',1);zoom on; grid on; hold on;
    H7 = plot(vecbinf,nineth_pow,'--k','Linewidth',1);zoom on; grid on; hold on;
    H8 = plot(vecbinf,th_pow,'--k','Linewidth',1);zoom on; grid on; hold on;
    %     H9 = plot(vecbinf,min_std,'--k','Linewidth',1);zoom on; grid on; hold on;
    %     H10 = plot(vecbinf,max_std,'--k','Linewidth',1);zoom on; grid on; hold on;
    % plot(vecbinf,HV1_b(:,2),'--b','Linewidth',1);zoom on; grid on; hold on;
    % plot(vecbinf,HV1_b(:,3),'--b','Linewidth',1);zoom on; grid on; hold on;
    % plot(vecbinf,HV1_c(:,2),'--m','Linewidth',1);zoom on; grid on; hold on;
    % plot(vecbinf,HV1_c(:,3),'--m','Linewidth',1);zoom on; grid on; hold on;

    %Axes log
    set(gca,'XScale', 'log', 'YScale', 'log');

    %Limites
    xlim([10^0 60]);
    ylim([0.2 60]);
    
    box on
    
    savename = "hv1_" + station + "." + location + "." + year; 
    saveas(gcf, fullfile(savedir, savename + ".jpg"))
    %  saveas(gcf,[savename '_norm.fig'])  

        
    %%%%%%%%%%%% TIME PLOT
    figure(2)
    imagesc(DVc,VECT_F,Stats_HV_mat)
    set(gca,'Yscale','log')
    set(gca,'Ydir','normal')
    % xlim([DVc(1) DVc(end)])
    ylim([fmin fmax])
    datetick('x','mm/yy','keeplimits')
    % xlabel('Date (mm/yy)')
    ylabel('Frequency (Hz)')
    load('lib/FFT_cmap.mat')
    colormap(cmap)
    c=colorbar;
    c.Label.String='HV';
    caxis([0 10])
    set(gca,'layer','top')
    box on
    % sdf('20_15')
    savename = "hv2_" + station + "." + location + "." + year;
    saveas(gcf, fullfile(savedir, savename + ".jpg"))
    %  saveas(gcf,[savename '_norm.fig'])  


    % in=datetime(2022, 5, 10);
    % in=datenum(in);
    % en=datetime(2022, 5, 20);
    % en=datenum(en);


    % HV normalizzato su max colonne
    Stats_HV_mat_norm=Stats_HV_mat./max(Stats_HV_mat);

    figure(3)
    imagesc(DVc,VECT_F,Stats_HV_mat_norm)
    % xlim([DVc(1) DVc(end)])
    set(gca,'Yscale','log')
    set(gca,'Ydir','normal')
    % xlim([in en])
    ylim([fmin fmax])
    datetick('x','mm/yy','keeplimits')
    % xlabel('Date (mm/yy)')
    ylabel('Frequency (Hz)')
    load('lib/FFT_cmap.mat')
    colormap(cmap);
    c=colorbar;
    c.Label.String='HV';
    caxis([0 1])
    set(gca,'layer','top')
    box on
    % sdf('20_15')

    savename = "hv3_" + station + "." + location + "." + year;
    saveas(gcf, fullfile(savedir, savename + ".jpg"))
    %  saveas(gcf,[savename '_norm.fig'])  
end