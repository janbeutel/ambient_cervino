

% Define the stations and their corresponding years
stations = {'MH44', 'MH36', 'MH38'};
years_by_station.MH36 = [2015 2016 2017 2018 2019 2020 2021 2022 2023 2024 2025];
years_by_station.MH38 = [2015 2016 2017 2018 2019];
years_by_station.MH44 = [2018 2019 2020 2021 2022 2023 2024 2025];

% Loop over each station and its years
for i = 1:length(stations)
    station = stations{i};
    years = years_by_station.(station);
    
    for j = 1:length(years)
        year = years(j);
        
        % Call the plotting functions
        CERVINO_plot_hilbert(year, station);
        CERVINO_plot_hv(year, station);
    end
end