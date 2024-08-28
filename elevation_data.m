%==========================================================================
% Dr. Muneer Al-Zubi
% Linkdin: https://www.linkedin.com/in/muneeralzubi85/
% Email:   muneermaz85@gmail.com

% For more details, the reader can refer to the following article:
% M M. Al-Zubi; M.-S. Alouini, "End-to-End Modelling and Simulation of NLOS Sub-6 GHz Backhaul via Diffraction for
% Internet Connectivity of Rural Areas", IEEE Open Journal of the Communications Society, 2023.
% https://ieeexplore.ieee.org/document/10306284
%==========================================================================

function [x z] = elevation_data(tx, rx, N_points, terrain_source)

%======input======
% tx              : tx site info
% rx              : rx site info
% N_points        : number of samples along path (resolution)
% terrain_source  : terrain elevation data source (e.g., SRTM 3 arc, 1 arc, etc )

% ======output======
% x                 : samples distance along path (m)
% z                 : elevation above sea level   (m)
%=================== 

d           = distance(tx,rx); % TX-RX distance
x_step      = d/N_points ; % distance step width (m)
x           = (0:N_points) * x_step ; % distance samples between TX & RX (m)

tx_loc      = [tx.Latitude tx.Longitude]; % TX location
rx_loc      = [rx.Latitude rx.Longitude]; % RX location

% Get elevation profile of terrain between tx and rx
[lat,lon]   = gcwaypts(tx_loc(1), tx_loc(2), rx_loc(1), rx_loc(2), N_points); % Get N_points lat/long points between TX and RX 
points      = rxsite("Name", "Test_Points", "Latitude", lat, "Longitude", lon,"AntennaHeight", zeros(1,N_points+1)); % create N_points sites along the path between TX & Rx
z           = elevation(points, 'Map', terrain_source)'; % get terrain elevation values at lat/long points
end