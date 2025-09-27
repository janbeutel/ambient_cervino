
% clc
close all 
clear all

year    = "2018";
network = "1I";
station = "MH44";
channel = "EHN.D";

data_directory = "../results/events/" + network + "/" + station + "/" + year + "/" + channel;
% List all .miniseed files in the directory
FileList1 = dir(fullfile(data_directory, '*.mat'));
TOT = size(FileList1,1);

DVec1=zeros(TOT,6);
% first=FileList(1).name;
% 
% year=str2num(first(7:10));
% month=str2num(first(12:13));
% day=str2num(first(15:16));
% hour=str2num(first(18:19));
% 
% DV_ref=[year month day hour];
% DV_ref(1,5:6)=0;



for ite= 1:1:TOT
    fprintf("Event %d, station %s, channel %s\n", ite, station, channel)

    % get the file name:
    folder = FileList1(ite).folder;
    filename = FileList1(ite).name;
    load(folder + "/" + filename); 

    L=length(SIG);

    %salta file danneggiat
    if L>1e6
        continue
    end

    mm=floor(sec/60);
    ss=floor(sec-(mm*60));

    year=str2double(filename(17:20));
    month=str2num(filename(21:22));
    day=str2num(filename(23:24));
    hour=str2num(filename(26:27));
    seconds=str2num(filename(33:36));
        
    DVec1(ite,:)=[year month day hour mm ss];
    
    SR=1/FS;   
 
    t=(X-X(1))*24*3600;
    sig=SIG*1.25;   %%correzione usando gain complessivo dei dati
    
    
    figure('color','w','Name',filename);

    subplot(4,1,1)
    plot(t,sig); 
    xlabel('Time [s]');
    ylabel ('V [m/s]');
    xlim([0 max(t)])
    hold on



    y=abs(fft(sig,L))/L;
    freq = FS*(0:(L/2))/L;
    sp(2:length(freq))=2*y(2:length(freq));      %true amplitude

    % freq = 0 : (1/(length(sig)*SR)) : 0.5/(SR);
    % sp = abs(fft(sig));
    % sp = sp(1:length(freq));

    subplot(4,1,2)
    plot (freq,sp)
    xlabel ('Frequency [Hz]');
    ylabel ('');
    xlim([0 125]);
    ylabel ('A [m/s]');
    hold on


    NF = fix(L/20); %piu' grande e' il numero piu' grandi gli intervalli in frequenza. numero di campioni in ogni finestra su cui calcolo fft
    OVERLP = fix (NF/1.01);
    % subplot(4,1,3:4);
    % hold on
    % specgram (sig,NF,FS,[],OVERLP);
    % xlabel ('Time (s)');
    % ylabel ('Frequency (Hz)');
    % xlim([0 max(t)])
    % ylim([0 125])
    % %colorbar('EastOutside')
    % box on
    % % colorbar('SouthOutside')
    % %caxis([-30 130])

    [a,f,tt]=specgram(sig,NF,FS,[],OVERLP);
    tt=tt+NF/2*SR;
    La=length(a(:,1));

    for i=1:length(a(1,:))
        a(:,i)=2*(abs(a(:,i))/La);
    end

    a=20*log10(a);

    subplot(4,1,3:4)
    hold on
    h = pcolor(tt,f,a);
    h.EdgeColor = 'none';
    %set(gca, 'YScale', 'log')
    ylim([0 125])
    %set(gca, 'YTick', [10 100 1000])
    xlim([t(1) t(end)])
    box on
    set(gca,'layer','top')
    xlabel('Time [s]')
    ylabel('Frequency [Hz]')
    colormap('jet')
    %c=colorbar
    %c.Label.String = 'Power [dB]';
    %caxis([-220 -100])
    %sdf('spectr')

    savedir = "../results/spectrograms/" + network + "/"+ station + "/" + year + "/" + channel;
    if ~exist(savedir, 'dir')  % check if folder exists
        mkdir(savedir);         % create folder if it doesn't exist
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    savename=[filename(1:36) + "_spectr"];
    % saveas(gcf,[savedir + "/" + filename + "_spectr.fig"])
    saveas(gcf,[savedir + "/" + savename + ".jpg"])
    clearvars -except DVec1 ite FileList1 year network channel station % SR FS
    close all
end


