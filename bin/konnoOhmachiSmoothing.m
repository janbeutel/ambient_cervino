function [smoothing_spectra]=konnoOhmachiSmoothing(spectra, frequencies, bandwidth, normalize);

% Create matrix to be filled with smoothing entries.
 smoothing_spectra = zeros(length(frequencies),1);
 for i=1:length(frequencies)
    freq=frequencies(i);
    window=konnoOhmachiSmoothingWindow(frequencies, freq, bandwidth, normalize);
    smoothing_spectra(i) = sum(window.*spectra);
 end
