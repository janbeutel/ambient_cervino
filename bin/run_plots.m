


% Define the stations and their corresponding years
NETWORK = '1I';
stations = {'MH36', 'MH38', 'MH44', 'MH54', 'MH48'};

years_by_station.MH36 = [2015 2016 2017 2018 2019 2020 2021 2022 2023 2024 2025];
years_by_station.MH38 = [2015 2016 2017 2018 2019];
years_by_station.MH44 = [2018 2019 2020 2021 2022 2023 2024 2025];
years_by_station.MH48 = [2019 2020 2021 2022 2023 2024 2025];
years_by_station.MH54 = [2019 2020];

% Define locations per station
locations_by_station.MH36 = {'A'};
locations_by_station.MH38 = {'A'};
locations_by_station.MH44 = {'B'};
locations_by_station.MH48 = {'A', 'B'};
locations_by_station.MH54 = {'A'};

% Loop over each station, year, and location
for i = 1:length(stations)
    station = stations{i};
    years = years_by_station.(station);
    locs = locations_by_station.(station);
    
    for j = 1:length(years)
        year = years(j);
        
        for k = 1:length(locs)
            location = locs{k};
            
            % Call the plotting functions
            CERVINO_plot_hilbert(NETWORK, year, station, location);
            CERVINO_plot_hv(NETWORK, year, station, location);
        end
    end
end