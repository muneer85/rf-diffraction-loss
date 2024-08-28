%==========================================================================
% Dr. Muneer Al-Zubi
% Linkdin: https://www.linkedin.com/in/muneeralzubi85/
% Email:   muneermaz85@gmail.com

% For more details, the reader can refer to the following article:
% M M. Al-Zubi; M.-S. Alouini, "End-to-End Modelling and Simulation of NLOS Sub-6 GHz Backhaul via Diffraction for
% Internet Connectivity of Rural Areas", IEEE Open Journal of the Communications Society, 2023.
% https://ieeexplore.ieee.org/document/10306284
%==========================================================================

function [Data1, Data2] = PathLoss (tx,rx,terrain_source, h_tree, N_points )
% Pathloss calaculation

% Free space loss (dB)
f         = tx.TransmitterFrequency * 1e-6; % MHz
d         = distance(tx,rx) * 1e-3; % distance tx-rx (km)
L_fs      = 32.5 + 20*log10(f) +20*log10(d); % d (km), f(Mhz)

% Pathloss using Knife-Edge (dB) Without radius
Dat1             = Single_Knife_Edge_Data(tx,rx,N_points,terrain_source, h_tree);
[L_KE theta]     = Single_Knife_Edge_Loss (tx,rx,Dat1);

% Pathloss correction term using rounded-obstacle with radius (dB)
F1        = Dat1(4);
Dat2     = Single_Rounded_Obs_Data(tx,rx,N_points,terrain_source, F1, h_tree);
Tmn       = Single_Rounded_Obs_Loss(tx, rx, Dat2);

% Total diffraction loss (dB) with obstacle radius
L_D       = L_KE + Tmn;

% Total loss (dB)
L         = L_D + L_fs;

Data1 = [L_fs, L_KE, L_D, L];
Data2 = [Dat2, d, theta];
%===========================================================

end