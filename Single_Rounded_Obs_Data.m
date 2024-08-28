%==========================================================================
% Dr. Muneer Al-Zubi
% Linkdin: https://www.linkedin.com/in/muneeralzubi85/
% Email:   muneermaz85@gmail.com

% For more details, the reader can refer to the following article:
% M M. Al-Zubi; M.-S. Alouini, "End-to-End Modelling and Simulation of NLOS Sub-6 GHz Backhaul via Diffraction for
% Internet Connectivity of Rural Areas", IEEE Open Journal of the Communications Society, 2023.
% https://ieeexplore.ieee.org/document/10306284
%==========================================================================

function G_info = Single_Rounded_Obs_Data(tx,rx,N_points,terrain_source, F1, h_tree)

%======input======
% tx              : tx site infor
% rx              : rx site info
% N_points        : number of samples along path (resolution)
% terrain_source  : terrain elevation data source (e.g., STRM 3 arc, 1 arc, etc )

% ======output======
% d1              :  distance from TX to intersection between (tx-tangent and rx-tangent)  
% d2              :  distance from RX to intersection between (tx-tangent and rx-tangent) 
% h               :  hight between vertex (intersection) and LOS-line (perpendicular) If the height is below this line, h is negative     
% R_curv          :  radius of curvature (m), e.q. 38 ITU-526
% F1              :  1st Fresnel zone radius at obstacle location. 
% F1_max          :  Max 1st Fresnel zone radius at the mid of path, d1 = d2

%=================== 

d           = distance(tx,rx); % TX-RX distance
x_step      = d/N_points ; % distance step width (m)
n_step      = N_points/d;  % 
f           = tx.TransmitterFrequency; % Center operating Freq. (Hz)
 
% Get elevation profile of terrain between TX and RX
[X Z] = elevation_data(tx, rx, N_points, terrain_source);

% location of TX/RX Antenna
h_tx        = tx.AntennaHeight; % tx antenna hight above ground
h_rx        = rx.AntennaHeight; % rx antenna hight above ground
TX_loc      = [X(1) Z(1)+h_tx]; % tx antenna location above sea level
RX_loc      = [X(N_points+1) Z(N_points+1)+h_rx]; % rx antenna location above sea level
x1          = TX_loc(1); %TX x
z1          = TX_loc(2); %TX z
x2          = RX_loc(1); %RX x
z2          = RX_loc(2); %RX z
%==============================
% Note: order of elements in the foolwing output arrays are reverse (from x-right to left)
% Get obstacles info (intersection points and lines for LOS, Fresnel Ellipse, obstacle top, etc)
[F1_max, Obs_Frex, Obs_Frez, Obs_losx, Obs_losz, Obs_top_loc, Obs_mid]= Get_Obstacle_Info(tx,rx,N_points,terrain_source, h_tree);
%==============================

n1               = length(Obs_Frex); 

if (n1>0) % if there is an intersection between terrain and (0.6 of 1st Fresnel zone)

% find the location of (top and mid-point @ LOS) of highest obstacle in the path 
[obs_max indx_obs1]            = max(Obs_top_loc(:,2));   % index of obstacle with highest top
obs_max_loc                    = Obs_top_loc(indx_obs1,:); % location of the highest top

[obs_max indx_obs2]            = max(Obs_mid(:,2));   % index of obstacle with highest top
Obs_max_mid                    = Obs_mid(indx_obs2,:); % mid point location (@ LOS line) of highest top

Obs_Frex_max                   = Obs_Frex(2*indx_obs1-1:2*indx_obs1); % Fresnel intersection x-coordinates with terrain for highest edge

% Obs_losx_max                   = Obs_losx(2*indx_obs1-1:2*indx_obs1); % LOS line intersection x-coordinates with terrain for highest edge

% Note: order of elements in the Obs_top_loc, Obs_mid, Obs_Frex, Obs_losx are reverse (from x-right to left)
%==============================

% Find radius of curvature

% method 1: ITU-526

% Rounded obstacle (fitting parabola to top of the obstacle)

% x              =  Obs_Frex_max(2):x_step:Obs_Frex_max(1); % distance samples between crossing point of Fresnel Ellipse and terrain
% z              =  Z(round(x*n_step));  % corresponding Z of these samples
% [vertex_x, vertex_z, tangent_tx_x, tangent_tx_z, tangent_rx_x, tangent_rx_z]  =  poly_fit(x,z,TX_loc,RX_loc); % fitting parabola to obstacle-top (see poly_fit function)

% %%====OR=========== more accurate====
x_Fre          =  Obs_Frex_max(2):x_step:Obs_Frex_max(1); % distance samples between crossing point of Fresnel Ellipse and terrain
z_Fre          =  Z(round(x_Fre*n_step));  % corresponding Z of these samples
z_test         =  (max(z_Fre) - z_Fre) < F1 ; % check if distance z_i<F1, see ITU-526 (rounded obstacle)
x              =  x_Fre(z_test);
z              =  z_Fre(z_test);
[vertex_x, vertex_z, tangent_tx_x, tangent_tx_z, tangent_rx_x, tangent_rx_z]  =  poly_fit(x,z,TX_loc,RX_loc); % fitting parabola to obstacle-top (see poly_fit function)

%=======================================

r_i     = 0;
X_n     = Obs_Frex_max(2):x_step:Obs_Frex_max(1) ; % distance samples between tx-rx (remove tx & rx location)
N       = 0;

for xx  = X_n  
zz      = Z(round(xx*n_step));
x_i     = obs_max_loc(1) - xx; % get xi, see ITU-526
z_i     = obs_max_loc(2) - zz; % get yi, see ITU-526

if (z_i <= F1) % check that max z_i is <= F1, see ITU-526
    if (z_i > 0) % avoid zero hight (e.g., at max elevation)
    r_i     = r_i + (x_i^2)/(2*z_i); % e.q. 38 ITU-526
    N       = N+1;
    end
end

end
R_curv  = r_i/N;     % radius of curvature, e.q. 38 ITU-526


% method 2: The occultation distance-based method






%==============================

% find intersection points between tx-tangent and rx-tangent lines (see ITU-526)
xx   =  0:x_step:x2;
m1  =  (z1-tangent_tx_z)/(x1-tangent_tx_x); % slop of tx-tangent line
m2  =  (z2-tangent_rx_z)/(x2-tangent_rx_x); % slop of rx-tangent line

Line_1       = m1 *(xx-x1) + z1; % equation of tx-tangent line
Line_2       = m2 *(xx-x2) + z2; % equation of rx-tangent line

[tangent_interc_x, tangent_interc_z]    = polyxpoly(xx, Line_1, xx, Line_2);  % intersection point between tx-tangent line and rx-tangent line
plot(tangent_interc_x,tangent_interc_z, 'x' , LineWidth=3, Color='g')         % plot this intersection point
%==============================

% find intersection point between (tangent_interc_x, tangent_interc_z) and LOS-line (perpendicular)

% First method (perpendicular)
s1              = (z2-z1)/(x2-x1); % slope of LOS line.
b               = z1-s1*x1; % LOS line between TX-RX: y=x.s1+b
s2              = -1/s1;% The slope of the perpendicular line is -1 over the original slope.
c               = tangent_interc_z - tangent_interc_x * s2 ; % h line from intersection-point to LOS-line : z=x.s2+c 
xIntersection   = (c - b) / (s1 - s2);  % x.s1+b = x.s2+c
zIntersection   = s2 * xIntersection + c;

% %Second method (direct)
% xIntersection    = tangent_interc_x;    % x coordinate of point @ LOS-line at cross-point between tx-tangent and rx-tangent 
% zIntersection    = s1*(xIntersection-x1)+z1,  % z coordinate of point @ LOS-line at obstacle-top location

%==============================

d1      = pdist2([x1 z1], [tangent_interc_x tangent_interc_z]); % distance from TX to intersection between (tx-tangent and rx-tangent) 
d2      = pdist2([x2 z2], [tangent_interc_x tangent_interc_z]); % distance from RX to intersection between (tx-tangent and rx-tangent) 
F1      = 17.3 * sqrt(1e-3* d1*d2/((f*1e-9)*(d1+d2))); % 1st Fresnel zone radius at highest obstacle location. d1,d2 are distances (km) from the TX/RX to obstruction. f in (GHz)
h       = pdist2([xIntersection zIntersection], [tangent_interc_x tangent_interc_z]); % hight between vertex (intersection) and LOS-line (perpendicular)
%h      = pdist2([Obs_max_mid(1) Obs_max_mid(2)], [tangent_interc_x tangent_interc_z]); % hight between top and LOS (mid point)

if (obs_max_loc(2) < zIntersection) % check if obstacle-top is under the LOS line
 h =  -h; % make h negative ==> theta is negative (angle of diffraction)
end
%======================================

hold on;
plot(Obs_mid(:,1),Obs_mid(:,2), 'x' , LineWidth=1, Color='r')       % plot mid-points of the top @ LOS line
plot(Obs_top_loc(:,1),Obs_top_loc(:,2),'*', LineWidth=1, Color='r') % plot point of obstacle-top 

plot([x1 tangent_interc_x], [z1 tangent_interc_z], LineWidth=1, Color='r') % plot line between TX and  vertex (Intersection)
plot([x2 tangent_interc_x], [z2 tangent_interc_z], LineWidth=1, Color='r') % plot line between RX and  vertex (Intersection)

plot([xIntersection tangent_interc_x], [zIntersection tangent_interc_z], LineWidth=1, Color='g') % plot h line between vertex (Intersection) and LOS

else

d1      = 0;
d2      = 0;
h       = 0;
F1      = 0;
R_curv  = 0;
end

G_info = [d1 d2 h R_curv F1 F1_max]; 
end
