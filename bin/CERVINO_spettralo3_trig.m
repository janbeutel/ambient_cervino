
% clc
close all 
clear all

year    = "2019";
network = "1I";
station = "MH44";
channel = "EHE.D";

data_directory = "../results/events/" + network + "/" + station + "/" + year + "/" + channel;
% List all .miniseed files in the directory
FileList1 = dir(fullfile(data_directory, '*.mat'));
TOT = size(FileList1,1);
DVec1=zeros(TOT,6);

savedir = "../results/spectrograms/" + network + "/"+ station + "/" + year + "/" + channel;
if ~exist(savedir, 'dir')  % check if folder exists
    mkdir(savedir);         % create folder if it doesn't exist
end

% Define log file path
logfile = "../results/spectrograms/" + network + "/" + station + "/" + year + "/" + channel + "_" + year + ".csv";
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
fprintf(fid, "eventfilename,detectiontime,windowlength,peak_amp,rms_amp,signalduration,year,month,day,hour,minute,second\n");


for ite= 1:1:TOT

    % get the file name:
    folder = FileList1(ite).folder;
    filename = FileList1(ite).name;
    load(folder + "/" + filename); 

    L=length(SIG);

    % Suppose tnum is a MATLAB datenum
    % tnum = 7.3847e+05;  
    dt = datetime(X(1), 'ConvertFrom', 'datenum');

    fprintf("Event %d, station %s, channel %s, %s\n", ite, station, channel, datestr(dt, 'yyyy-mm-dd HH:MM:SS'))

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
    
    % signal duration
    thr = 0.1 * max(abs(sig));   % 10% of max amplitude
    idx = find(abs(sig) > thr);
    duration = t(idx(end)) - t(idx(1));
    
    y=abs(fft(sig,L))/L;
    freq = FS*(0:(L/2))/L;
    sp(2:length(freq))=2*y(2:length(freq));      %true amplitude

    NF = fix(L/20); %piu' grande e' il numero piu' grandi gli intervalli in frequenza. numero di campioni in ogni finestra su cui calcolo fft
    OVERLP = fix (NF/1.01);

    [a,f,tt]=specgram(sig,NF,FS,[],OVERLP);
    tt=tt+NF/2*SR;
    La=length(a(:,1));

    for i=1:length(a(1,:))
        a(:,i)=2*(abs(a(:,i))/La);
    end

    a=20*log10(a);


    figure('color','w','Name',filename);
    h1 = subplot(4,1,2);
    plot(t,sig); 
    xlabel('Time [s]');
    ylabel ('V [m/s]');
    xlim([0 max(t)])
    hold on
    set(h1, 'Position', [0.1 0.5 0.85 0.15]);  % [left bottom width height]

    h2 = subplot(4,1,1);
    plot (freq,sp)
    xlabel ('Frequency [Hz]');
    ylabel ('');
    xlim([0 125]);
    ylabel ('A [m/s]');
    hold on
    set(h2, 'Position', [0.1 0.75 0.85 0.15]);  % [left bottom width height]

    h3 = subplot(4,1,3:4);
    h = pcolor(tt,f,a);
    h.EdgeColor = 'none';
    ylim([0 125])
    xlim([t(1) t(end)])
    box on
    set(gca,'layer','top')
    xlabel('Time [s]')
    ylabel('Frequency [Hz]')
    colormap('jet')
    hold on
    set(h3, 'Position', [0.1 0.1 0.85 0.3]);  % [left bottom width height]

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    savename=[filename(1:36) + "_spectr"];

    % Write line to logfile
    fprintf(fid, "%s,%s,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n", filename(1:36), datestr(dt, 'yyyy-mm-dd HH:MM:SS'), L/FS, max(abs(sig)), rms(sig), duration, year, month, day, hour, mm, ss);

    % saveas(gcf,[savedir + "/" + filename + "_spectr.fig"])
    saveas(gcf,[savedir + "/" + savename + ".jpg"])
    clearvars -except DVec1 ite FileList1 year network channel station fid savedir % SR FS
    close all
end

% Close logfile
fclose(fid);
