function [smoothing_window]=konnoOhmachiSmoothingWindow(frequencies,center_frequency,bandwidth,normalize)

% If the center_frequency is 0 return an array with zero everywhere except at zero.
if center_frequency == 0
    smoothing_window = zeros(1,length(frequencies));
    smoothing_window(frequencies == 0.0)=1.0;
else
    % Calculate the bandwidth*log10(f/f_c)
    smoothing_window=bandwidth*log10(frequencies/center_frequency);
    % Just the Konno-Ohmachi formulae.
    smoothing_window=(sin(smoothing_window)./smoothing_window).^4;
    % Check if the center frequency is exactly part of the provided
    % frequencies. This will result in a division by 0. The limit of f->f_c is
    % one.
    smoothing_window(frequencies == center_frequency) = 1.0;
    % Also a frequency of zero will result in a logarithm of -inf. The limit of
    % f->0 with f_c!=0 is zero.
    smoothing_window(frequencies == 0.0) = 0.0;
    % Normalize to one if wished.
    if normalize
        smoothing_window=smoothing_window./sum(smoothing_window);
    end
end
