%==========================================================================
% Dr. Muneer Al-Zubi
% Linkdin: https://www.linkedin.com/in/muneeralzubi85/
% Email:   muneermaz85@gmail.com

% For more details, the reader can refer to the following article:
% M M. Al-Zubi; M.-S. Alouini, "End-to-End Modelling and Simulation of NLOS Sub-6 GHz Backhaul via Diffraction for
% Internet Connectivity of Rural Areas", IEEE Open Journal of the Communications Society, 2023.
% https://ieeexplore.ieee.org/document/10306284
%==========================================================================

function Tmn = Single_Rounded_Obs_Loss(tx, rx, G_info)

% see ITU-526

%======input======
% tx              : tx site infor
% rx              : rx site info
% G_info          : values of d1, d2, h, R

% ======output======
% Tmn             : correction pathloss due to rounded-obstacle with radius (dB)
%=================== 

f       = tx.TransmitterFrequency; % Center operating Freq. (Hz)
C       = physconst('light');      % Light Speed in vacuum (m/s)
Lambda  = C/f;                     % Wavelength (m)
d1      = G_info(1);               % distances of the tx from the vertex (intersection point)             
d2      = G_info(2);               % distances of the rx from the vertex (intersection point) 
h       = G_info(3);               % hight between vertex (intersection) and LOS-line (perpendicular). If the height is below this line, h is negative     
R       = G_info(4);               % Obstacle radius of curvature (m)

m       = R * ((d1+d2)/(d1*d2)) / ((pi*R/Lambda)^(1/3)); 
n       = h * ((pi*R/Lambda)^(2/3))/R;

test    = m*n;

if (test <= 4)
Tmn     = 7.2 * m^0.5 - (2-12.5*n)*m  + 3.6*m^(3/2) - 0.8*m^2;  % e.q. 34a ITU 526
else
Tmn     =  -6 - 20*log10(m*n) + 7.2*m^0.5 - (2-17*n)*m + 3.6*m^(3/2) - 0.8*m^2;  % e.q. 34b ITU 526
end  

if (h==0)
  Tmn  =0;
end
end