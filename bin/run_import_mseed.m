function run_import_mseed(year, station, channel)
    fprintf("Running detection for year %d, station %s, channel %s\n", year, station, channel);
    CERVINO_import_mseed_orari(year, station, channel);
end
