network  = "1I";
year     = "2019";
station  = "MH48";
location = "B";
channel  = "EHE.D";

fprintf("Running detection for network %s, year %s, station %s, location %s, channel %s\n", network, year, station, location, channel);
CERVINO_import_mseed_orari(network, year, station, location, channel);
% CERVINO_import_mseed(year, station, channel);
