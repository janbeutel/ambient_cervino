year    = "2018"
station = "MH44"
channel = "EHE.D"

fprintf("Running detection for year %s, station %s, channel %s\n", year, station, channel);
CERVINO_import_mseed_orari(year, station, channel);
% CERVINO_import_mseed(year, station, channel);
