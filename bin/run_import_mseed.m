network  = "1I";
year     = "2015";
station  = "MH36";
location = "A";
channel  = "EHE.D";

fprintf("Running detection for network %s, year %s, station %s, location %s, channel %s\n", network, year, station, location, channel);
CERVINO_import_mseed_orari(network, year, station, location, channel);
% CERVINO_import_mseed(year, station, channel);
