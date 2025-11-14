

% Define the stations and their corresponding years
NETWORK=1I
stations = {'MH54', 'MH48', 'MH44', 'MH36', 'MH38'};
years_by_station.MH36 = [2015 2016 2017 2018 2019 2020 2021 2022 2023 2024 2025];
years_by_station.MH38 = [2015 2016 2017 2018 2019];
years_by_station.MH44 = [2018 2019 2020 2021 2022 2023 2024 2025];
years_by_station.MH48 = [2019 2020 2021 2022 2023 2024 2025];
years_by_station.MH54 = [2019 2020];

% Loop over each station and its years
for i = 1:length(stations)
    station = stations{i};
    years = years_by_station.(station);
    
    for j = 1:length(years)
        year = years(j);
        
        % Call the plotting functions
        CERVINO_plot_hilbert(network, year, station, location);
        CERVINO_plot_hv(network, year, station, location);
    end
end