%==========================================================================
% Dr. Muneer Al-Zubi
% Linkdin: https://www.linkedin.com/in/muneeralzubi85/
% Email:   muneermaz85@gmail.com

% For more details, the reader can refer to the following article:
% M M. Al-Zubi; M.-S. Alouini, "End-to-End Modelling and Simulation of NLOS Sub-6 GHz Backhaul via Diffraction for
% Internet Connectivity of Rural Areas", IEEE Open Journal of the Communications Society, 2023.
% https://ieeexplore.ieee.org/document/10306284
%==========================================================================

function [L_KE theta] = Single_Knife_Edge_Loss(tx, rx, G_info)

% Fresnel-Kirchoff loss

%======input======
% tx              : tx site infor
% rx              : rx site info
% G_info          : values of d1, d2, h 

% ======output======
% L_KE                 : pathloss due to knife-edge obstacle (dB)
%=================== 
d           = distance(tx,rx); % TX-RX distance
f           = tx.TransmitterFrequency; % Center operating Freq. (Hz)
C           = physconst('light'); % Light Speed in vacuum (m/s)
Lambda      = C/f; % Wavelength (m)
d1          = G_info(1); % distances between tx and obstacle top               
d2          = G_info(2); % distances between rx and obstacle top
h           = G_info(3); % height of obstacle top above the straight line joining tx & rx.
                         % If the height is below this line, h is negative          
alph1       = G_info(6); % distances between rx and obstacle top
alph2       = G_info(7); % distances between rx and obstacle top


%========== alternative method to find v =========
%calculate v from angles alpha1 and alpha2                                            % other end. it has  sign of h.
%v=sqrt((2*d/Lambda)*alph1*alph2);
%===================

v           = h*sqrt(2/Lambda*((1/d1)+(1/d2))); % calculate v from h, d1, d2
theta       = (v^2)*Lambda/(2*h) * 180/pi; % diffraction angle (degree). its sign is the same as that of h. The angle is assumed
                                           % to be less than about 0.2 rad, or roughly 12º

% real and imaginary part of the Complex Fresnel Integral (eqs. 7, 30, ITU-526)
C_fun       = @(s) cos((pi*s.^2)/2);
S_fun       = @(s) sin((pi*s.^2)/2);
C           = integral(C_fun,0,v);
S           = integral(S_fun,0,v);
Fun1        = sqrt((1-C-S)^2+(C-S)^2); 

L_KE        = -20*log10(Fun1/2); % Diffracion loss (dB)

if (h==0)
  L_KE  =0;
end

end