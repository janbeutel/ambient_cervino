close all
clear all
clc

year    = "2019";
network = "1I";
station = "MH44";

data_directory = "../results/fft/" + network + "/" + station + "/" + year + "/";
% List all .miniseed files in the directory
FileList = dir(fullfile(data_directory, '*.mat'));
TOT = size(FileList,1);

savedir = "../results/plots/";
if ~exist(savedir, 'dir')  % check if folder exists
    mkdir(savedir);         % create folder if it doesn't exist
end

fprintf("Running CERVINO_plot_hibert for year %s, station %s, %d\n", year, station, TOT)


channel = "EHE.D";
filename=["fft_" + station + "_" + channel + "_" + year + ".mat"]
load(data_directory + filename) 
FFT_E=FFT_all;
DV_E=DVec;

channel = "EHN.D";
filename=["fft_" + station + "_" + channel + "_" + year + ".mat"];
load(data_directory + filename) 
FFT_N=FFT_all;
DV_N=DVec;

channel = "EHZ.D";
filename=["fft_" + station + "_" + channel + "_" + year + ".mat"];
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



channel = "EHN.D";

data_directory = "../results/classification/" + network + "/" + station + "/" + year + "/";
% List all .miniseed files in the directory
FileList = dir(fullfile(data_directory, '*.mat'));
TOT = size(FileList,1);

filename=["classification_" + station + "_" + channel + "_" + year + ".mat"]
load(data_directory + filename) 

% cd('C:\Users\Valer\OneDrive - Politecnico di Torino\CERVINO RISULTATI\2015-2025\MH36\microseismicity')
% load('CLASS_MH44_N.mat');   %%%%%%%%%%%%%%%
load('fieldwork_dates.mat');
% load('rockfall_dates.mat');
% load('helicopter_dates.mat');
% load('MH27_temperature.mat');
% load('MH11_temperature.mat'); 
% time(1)=[];

data=Hilbert;
Dv=datenum(DVec1);


start_fw=datenum(start_fw);
end_fw=datenum(end_fw);

% start_rf=datenum(start_rf);
% end_rf=datenum(end_rf);

% start_eli=datenum(start_eli);
% end_eli=datenum(end_eli);


minDv=min(Dv);
maxDv=max(Dv);

% time=datenum(time);
% idx=(time >= minDv & time <= maxDv);
% time=time(idx);
% T_100cm=T_100cm(idx);

% idx2=(DVc >=minDv & DVc<=maxDv);
% DVc=DVc(idx);
% Stats_HV_mat=Stats_HV_mat(:,idx);

% DVec=DVec1;
% a=find(Hilbert(:,9)>35);
% data(a,:)=[];
% Dv=datenum(DVec);
% Dv(a)=[];


data_f=data;
Dv_f=Dv;
a=find(Hilbert(:,9)>35);
data_f(a,:)=[];
Dv_f(a)=[];


figure (1)
hold on 
% subplot(4,1,1)
% plot(time, T_100cm)
% xlim([minDv maxDv])
% ylim([min(T_100cm) max(T_100cm)])
% xlabel('Date (mm/yy)')
% ylabel('T (°C)')
% datetick('x','mm/yy', 'keeplimits')

subplot(3,1,1)
plot(Dv, data(:,9),'.')
datetick('x','mm/yy', 'keeplimits')
xlim([minDv maxDv])
xlabel('Date (mm/yy)')
ylabel('F(Hz)')
hold on
for i = 1:length(start_fw)
    idx = Dv >= start_fw(i) & Dv <= end_fw(i);
    plot(Dv(idx), data(idx,9), 'r.', 'MarkerSize', 12)
end

subplot(3,1,2)
plot(Dv_f, data_f(:,9),'.')
datetick('x','mm/yy', 'keeplimits')
xlim([minDv maxDv])
xlabel('Date (mm/yy)')
ylabel('F(Hz)')
hold on
hold off


subplot(3,1,3)
imagesc(DVc,VECT_F,Stats_HV_mat)
set(gca,'Yscale','log')
set(gca,'Ydir','normal')
xlim([minDv maxDv])
ylim([5 50])
datetick('x','mm/yy','keeplimits')
xlabel('Date (mm/yy)')
ylabel('Frequency (Hz)')
load('FFT_cmap.mat')
colormap(cmap);
c=colorbar;
c.Label.String='HV'
caxis([0 10])
set(gca,'layer','top')
box on
% sdf('20_15')




figure(2)
hold on
% subplot(5,1,1)
% plot(time, T_100cm)
% xlim([minDv maxDv])
% ylim([min(T_100cm) max(T_100cm)])
% xlabel('Date (mm/yy)')
% ylabel('T (°C)')
% datetick('x','mm/yy', 'keeplimits')
% 
% subplot(5,1,2)
% plot(Dv, data(:,2),'.')
% datetick('x','mm/yy', 'keeplimits')
% xlim([minDv maxDv])
% xlabel('Date (mm/yy)')
% ylabel('Kurtosis (-)')
% hold on
% for i = 1:length(start_fw)
%     idx = Dv >= start_fw(i) & Dv <= end_fw(i);
%     plot(Dv(idx), data(idx,2), 'r.', 'MarkerSize', 12)
% end

subplot(3,1,1)
plot(Dv, data(:,4),'.')
datetick('x','mm/yy', 'keeplimits')
xlim([minDv maxDv])
xlabel('Date (mm/yy)')
ylabel('Uniform duration (s)')
hold on
for i = 1:length(start_fw)
    idx = Dv >= start_fw(i) & Dv <= end_fw(i);
    plot(Dv(idx), data(idx,4), 'r.', 'MarkerSize', 12)
end

subplot(3,1,2)
plot(Dv, data(:,10),'.')
datetick('x','mm/yy', 'keeplimits')
xlim([minDv maxDv])
ylim([0,0.0002])
xlabel('Date (mm/yy)')
ylabel('Max amplitude (dB)')
hold on
for i = 1:length(start_fw)
    idx = Dv >= start_fw(i) & Dv <= end_fw(i);
    plot(Dv(idx), data(idx,10), 'r.', 'MarkerSize', 12)
end

subplot(3,1,3)
plot(Dv, data(:,1),'.')
datetick('x','mm/yy', 'keeplimits')
xlim([minDv maxDv])
xlabel('Date (mm/yy)')
ylabel('Amax/Amean')
hold on
for i = 1:length(start_fw)
    idx = Dv >= start_fw(i) & Dv <= end_fw(i);
    plot(Dv(idx), data(idx,1), 'r.', 'MarkerSize', 12)
end
% sdf('20_15')