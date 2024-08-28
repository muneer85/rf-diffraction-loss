%==========================================================================
% Dr. Muneer Al-Zubi
% Linkdin: https://www.linkedin.com/in/muneeralzubi85/
% Email:   muneermaz85@gmail.com

% For more details, the reader can refer to the following article:
% M M. Al-Zubi; M.-S. Alouini, "End-to-End Modelling and Simulation of NLOS Sub-6 GHz Backhaul via Diffraction for
% Internet Connectivity of Rural Areas", IEEE Open Journal of the Communications Society, 2023.
% https://ieeexplore.ieee.org/document/10306284
%==========================================================================

function [F1_max, Obs_Frex, Obs_Frez, Obs_losx, Obs_losz, Obs_top_loc, Obs_mid] = Get_Obstacle_Info(tx,rx,N_points,terrain_source, h_tree)

%======input======
% tx              : tx site infor
% rx              : rx site info
% N_points        : number of samples along path (resolution)
% terrain_source  : terrain elevation data source (e.g., STRM 3 arc, 1 arc, etc )

% ======output======
% F1_max                : Max 1st Fresnel zone radius (m) @ mid of the path 
% [Obs_Frex, Obs_Frez]  : intersection point of 0.6 of Fresnel Ellipse and terran elevation 
% Obs_top_loc           : locations of all knife-edge top points
% Obs_mid               : locations of all mid-points @ LOS line for knife-edge tops
%=================== 

d           = distance(tx,rx); % TX-RX distance
x_step      = d/N_points ;     % distance step width (m)
n_step      = N_points/d;
f           = tx.TransmitterFrequency; % Center operating Freq. (Hz)

% Get elevation profile of terrain between TX and RX
[X Z] = elevation_data(tx, rx, N_points, terrain_source);

% location of TX/RX Antenna
h_tx        = tx.AntennaHeight; % tx antenna hight above ground
h_rx        = rx.AntennaHeight; % rx antenna hight above ground
TX_loc      = [X(1) Z(1)+h_tx]; % tx antenna location above sea level
RX_loc      = [X(N_points+1) Z(N_points+1)+h_rx]; % rx antenna location above sea level
x1          = TX_loc(1); % TX x
z1          = TX_loc(2); % TX z
x2          = RX_loc(1); % RX x
z2          = RX_loc(2); % RX z
%==============================

% Fresnel Ellipse 
F1_max  = 8.657 * sqrt(d*1e-3/(f*1e-9)); % Max 1st Fresnel zone radius (m) @ mid of the path 

a       = 1/2*sqrt((x2-x1)^2+(z2-z1)^2); % (major_axis)/2
b       = F1_max * 0.60; % (second_axis)/2 [60% of Fresnel zoon should be cleare]

t       = linspace(0,-1*pi,N_points+1); % 0 to -pi give the lower half of the Ellipse
xa      = a*cos(t);
zb      = b*sin(t); 
w       = atan2(z2-z1,x2-x1); % angle of two points
Xe      = (x1+x2)/2 + xa*cos(w) - zb*sin(w); % Fresnel Ellipse X
Ze      = (z1+z2)/2 + xa*sin(w) + zb*cos(w); % Fresnel Ellipse Z
%==============================

% Find intersection points between terrain & LOS   and    between Terrain & Fresnel Ellipse
Line_LOS            = linspace(TX_loc(2),RX_loc(2),N_points+1); % LOS line TX-RX

[Obs_Frex Obs_Frez] = polyxpoly(X, Z, Xe, Ze);      % intersection point of 0.6 of Fresnel Ellipse and terran elevation 
[Obs_losx Obs_losz] = polyxpoly(X, Z, X, Line_LOS); % intersection point of LOS line and terran elevation
%==============================

% Find top and (mid-point @ LOS) of each knife edge

n1               = length(Obs_Frex); 

% Find top of each knife edge
Obs_top_loc     = zeros(n1/2, 2);  % each two cross-points give one top-point       
j1=0;
for i=1:2:n1
j1= j1+1;
x_obs                         = [Obs_Frex(i+1):x_step:Obs_Frex(i)]; % add +1 to a void 0 value
[Z_obs  indx_obs]             = max(Z(round(x_obs*n_step))); 
Obs_top_loc(j1,:)             = [x_obs(indx_obs) Z_obs+h_tree]; 
%Obs_mid(j1,:)                = [sum(Obs_Frex(i:(i+1)))/2 sum(Obs_Frey(i:(i+1)))/2];
end

% Find mid-point @ LOS for each knife edge
n1= length(Obs_losx); 
Obs_mid         = zeros(n1/2, 2); % each two cross-points give one mid-point
j1              = 0;
for i=1:2:n1
j1              = j1+1;
x_obs           = [Obs_losx(i+1):x_step:Obs_losx(i)];
Obs_mid(j1,:)   = [sum(Obs_losx(i:(i+1)))/2 sum(Obs_losz(i:(i+1)))/2];
end


if (n1>0)   % if there is an intersection between terrain and 0.6 of 1st Fresnel zone

plot(Obs_losx, Obs_losz, 'x' , LineWidth=1, Color='r') % plot cross-points between LOS and terrain
plot(Obs_Frex, Obs_Frez, 'x' , LineWidth=1, Color='r') % plot cross-points between 0.6 of Fresnel Ellipse and terrain

else

disp("Clear Fresnel Zoon") % if no cross between terrain and 0.6 of Fresnel Ellipse

end

hold on;
plot(X, Z, LineWidth=2, Color='k') % plot elevation
plot(X, Line_LOS, LineWidth=2, Color='b') % plot LOS line between TX-RX
plot(Xe, Ze, LineWidth=2, Color='y') % plot lower-half of Fresnel Ellipse
plot([x1 x1], [Z(1) z1], LineWidth=2, Color='r', Marker='_') % plot TX
plot([x2 x2], [Z(N_points+1) z2], LineWidth=2, Color='r', Marker='_') % plot RX

xlim([0 d])
title('Elevation Between TX and RX');
xlabel('Distance (m)') 
ylabel('Elevation (m)')

% axis equal %show perpendicular angle between h and LOS

% Plot points along the LOS path in SiteViewer
% Values   = ones(length(lat_test),1);
% Values    = z;
% tbl       = table(lat, lon, Values);
% Data      = propagationData(tbl);
% plot(Data);

end