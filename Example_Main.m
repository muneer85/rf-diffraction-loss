%==========================================================================
% Dr. Muneer Al-Zubi
% Linkdin: https://www.linkedin.com/in/muneeralzubi85/
% Email:   muneermaz85@gmail.com

% For more details, the reader can refer to the following article.
% If you use this code in your work, please cite the following article.

% M. M. Al-Zubi and M.-S. Alouini, "End-to-End Modelling and Simulation of NLOS Sub-6 GHz Backhaul via
% Diffraction for Internet Connectivity of Rural Areas," IEEE Open Journal of the Communications Society, 2023.
% https://ieeexplore.ieee.org/document/10306284
%==========================================================================

clear
clc
close all force

%==========================================================================
% define the terrain elevation data source (i.e., digital terrian/elevation model)
terrain_name          = 'gmted2010'; % default digital terrain model is (GMTED2010). You can replace it by any other names if you will use 
                                     % an external Digital Terrain Elevation Data (DTED) files of other models such as SRTM 
                                     % (i.e., *.dt0, *.dt1, or *.dt2) using the fellowing command: 

dtedfiles             = ["terrain_files/file1.dt2"];   % name and path of the terrain source file, where multiple files can                                                  
                                                       % be added, e.g.,  ["file1" "file2" ... "fileN"]                                             
try 
    addCustomTerrain(terrain_name, dtedfiles, "FillMissing", true) % change the terrain data source
catch ME 
    if strcmp(ME.identifier, 'shared_terrain:terrain:TerrainNameExists') 
        disp("The terrain file is already uploaded.")
        %Do nothing 
    else 
        rethrow(ME); 
    end 
end 
%==========================================================================
% define the location of transmitter and receiver
TX_Lat                       = [40.375111111111110]; % TX site latitude 
TX_Lon                       = [-4.258472222222222]; % TX site longitude
RX_Lat                       = [40.375833333333330]; % RX site latitude 
RX_Lon                       = [-4.396055555555556]; % RX site longitude

viewer                       = siteviewer("Terrain",terrain_name); % create site viewer
       
tx                           = txsite("Name", "TX", "Latitude", TX_Lat, "Longitude", TX_Lon); % create TX site
tx.TransmitterFrequency      = 5800e6; % set the transmitter frequency in Hz
tx.AntennaHeight             = 30;     % set the transmitter (e.g., BS) antenna height in m

rx                           = rxsite("Name", "RX", "Latitude", RX_Lat, "Longitude", RX_Lon); % create RX site
rx.AntennaHeight             = 2; % set the receiver (e.g., mobile) antenna height in m

h_tree                       = 0; % set an average value for the tree height in the path between TX and RX
N_points                     = 5000; % set the number of points used for sampling the path between TX and RX. 
                                     % the larger number provides higher resolution but with more computational time
show(tx) % show TX site in the site viewer
show(rx) % show RX site in the site viewer
los(tx,rx) % show the LOS/NLOS link between TX and RX
%==========================================================================

% Pathloss Calculation
[data1, data2]        = PathLoss (tx,rx,terrain_name, h_tree, N_points);

disp("L_fs     = " + data1(1) + " dB")  % L_fs:   free space loss
disp("L_KE     = " + data1(2) + " dB")  % L_KE:   Knife edge diffraction loss
disp("L_D      = " + data1(3) + " dB")  % L_D:    rounded obstacle diffraction loss
disp("L        = " + data1(4) + " dB")  % L:      total pathloss [L_fs + L_D]

disp("d1       = " + data2(1) + " m")   % d1:     the distances from the TX antenna to the top of the obstacle
disp("d2       = " + data2(2) + " m")   % d2:     the distances from the RX antenna to the top of the obstacle
disp("h        = " + data2(3) + " m")   % h:      the distance from the top of the obstacle to the LOS line that connects the two antennas
                                        %         (h is negative if the obstacle is below this line)
disp("R        = " + data2(4) + " m")   % R:      the obstacle radius of curvature in m
disp("F1       = " + data2(5) + " m")   % F1:     1st Fresnel zone radius at obstacle location.
disp("F1_max   = " + data2(6) + " m")   % F1_max: maximum 1st Fresnel zone radius at the mid of path, d1 = d2
disp("d        = " + data2(7) + " km")  % d:      distance between TX and RX in Km
disp("theta    = " + data2(8) + " deg") % theta:  diffraction angle in degree
