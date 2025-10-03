%Fonction: func_spectre_puissance
%Auteur: Pierre Bottelin, sept 2011

%Description: Calcule a partir d'un signal temporel, le spectre de puissance du signal
%Avec normalisation par la longueur de fenetre temporelle
%Avec redistribution de l'énergie totale sur [0 fs/2]

%Rqs: 
%- nfft courant: nfft=2^nextpow2(length(signal))
%puissance de 2 pour optimisation algo FFT matlab
%- apodisation courante: 0.1 (=10%) cosine taper (=tukeywin)

function [Vect_f,Spectre_puissance]=func_spectre_puissance(signal,fs,nfft,apodisation);

%Mise a la bonne dimension
if size(signal,1)>1
    signal=signal';
end

%-----------APODISATION------------------------------------
signal_apod=signal.*(tukeywin(length(signal),apodisation)');

%-----------VECTEUR FREQUENCE-----------------------------
Vect_f=(fs/2)*linspace(0,1,nfft/2+1);

%-----------SPECTRE----------------------------------------
Spectre_complexe=fft(signal_apod,nfft)/length(signal);
% Spectre_puissance=2*abs(Spectre_complexe(1:nfft/2+1));
Spectre_puissance=abs(Spectre_complexe(1:nfft/2+1));

clear signal fs apodisation signal_apod nfft Spectre_complexe


