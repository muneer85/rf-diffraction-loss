%==========================================================================
% Dr. Muneer Al-Zubi
% Linkdin: https://www.linkedin.com/in/muneeralzubi85/
% Email:   muneermaz85@gmail.com

% For more details, the reader can refer to the following article:
% M M. Al-Zubi; M.-S. Alouini, "End-to-End Modelling and Simulation of NLOS Sub-6 GHz Backhaul via Diffraction for
% Internet Connectivity of Rural Areas", IEEE Open Journal of the Communications Society, 2023.
% https://ieeexplore.ieee.org/document/10306284
%==========================================================================

function  G_info = Single_Knife_Edge_Data(tx,rx,N_points,terrain_source, h_tree)

%======input======
% tx              : tx site infor
% rx              : rx site info
% N_points        : number of samples along path (resolution)
% terrain_source  : terrain elevation data source (e.g., STRM 3 arc, 1 arc, etc )

% ======output======
% d1              :  distances between tx and obstacle top 
% d2              :  distances between rx and obstacle top 
% h               :  height of obstacle top above the straight line joining tx & rx. If the height is below this line, h is negative     
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

% Get obstacles info (intersection points and lines for LOS, Fresnel Ellipse, obstacle top, etc)
[F1_max, Obs_Frex, Obs_Frez, Obs_losx, Obs_losz, Obs_top_loc, Obs_mid]= Get_Obstacle_Info(tx,rx,N_points,terrain_source, h_tree);

%==============================
n1               = length(Obs_Frex); 

if (n1>0) % if there is an intersection between terrain and 0.6 of 1st Fresnel zone

% find the location of (top and mid-point @ LOS) of highest obstacle in the path
[obs_max indx_obs]            = max(Obs_top_loc(:,2));   % index of obstacle with highest top
obs_max_loc                   = Obs_top_loc(indx_obs,:); % location of the highest top

[obs_max indx_obs]            = max(Obs_mid(:,2));   % index of obstacle with highest top
Obs_max_mid                   = Obs_mid(indx_obs,:); % mid point location (@ LOS line) of highest top
%==============================

% Find intersection point between max-top and LOS-line @ perpendicular 

% first method (perpendicular)
s1               = (z2-z1)/(x2-x1); % Determine slope of LOS line
b                = z1-s1*x1;        % LOS line between TX-RX: z=x.s1+b
s2               = -1/s1;           % slope of the perpendicular line is -1 over the original slope.
c                = obs_max_loc(2) - obs_max_loc(1) * s2 ; % h line from obstacle-top to LOS line : z=x.s2+c 
xIntersection    = (c - b) / (s1 - s2);  % x.s1+b = x.s2+c
zIntersection    = s2 * xIntersection + c; 

% % second method (direct)
% xIntersection    = obs_max_loc(1);    % x coordinate of point @ LOS-line at obstacle-top location
% zIntersection    = s1*(xIntersection-x1)+z1,  % z coordinate of point @ LOS-line at obstacle-top location

%==============================

d1      = pdist2([x1 z1], [obs_max_loc(1) obs_max_loc(2)]); % distance from TX to obstacle-top
d2      = pdist2([x2 z2], [obs_max_loc(1) obs_max_loc(2)]); % distance from RX to obstacle-top
F1      = 17.3 * sqrt(1e-3* d1*d2/((f*1e-9)*(d1+d2))); % 1st Fresnel zone radius at highest obstacle location. d1,d2 are distances (km) from the TX/RX to obstruction. f in (GHz)

h       = pdist2([xIntersection zIntersection], [obs_max_loc(1) obs_max_loc(2)]);  % hight between obstacle-top and LOS (perpendicular)
% h       = pdist2([Obs_max_mid(1) Obs_max_mid(2)], [obs_max_loc(1) obs_max_loc(2)]); % hight between obstacle-top and LOS (mid point)
if (obs_max_loc(2) < zIntersection) % check if obstacle-top is under the LOS line
 h =  -h; % make h negative ==> theta is negative (angle of diffraction)
end
%==============================

% find angles alpha1 and alpha2
x3=obs_max_loc(1);
z3=obs_max_loc(2);
v_1 = [x2,z2,0] - [x1,z1,0];
v_2 = [x3,z3,0] - [x1,z1,0];
alph1 = atan2(norm(cross(v_1, v_2)), dot(v_1, v_2)); % angles in radians between the top of the obstacle and tx as seen from the
                                                     %  other end. it has  sign of h.
v_1 = [x1,z1,0] - [x2,z2,0];
v_2 = [x3,z3,0] - [x2,z2,0];
alph2 = atan2(norm(cross(v_1, v_2)), dot(v_1, v_2)); % angles in radians between the top of the obstacle and rx as seen from the
         
%======================================

hold on;
plot(Obs_mid(:,1),Obs_mid(:,2), 'x' , LineWidth=1, Color='r') % plot mid-points of the top @ LOS line
plot(Obs_top_loc(:,1),Obs_top_loc(:,2),'*', LineWidth=1, Color='r') % plot point of obstacle-top 

plot([x1 obs_max_loc(1)], [z1 obs_max_loc(2)], LineWidth=1, Color='r') % plot line between TX and  obstacle-top
plot([x2 obs_max_loc(1)], [z2 obs_max_loc(2)], LineWidth=1, Color='r') % plot line between RX and  obstacle-top

plot([xIntersection obs_max_loc(1)], [zIntersection obs_max_loc(2)], LineWidth=1, Color='g')  % plot h line between obstacle-top and LOS (perpendicular)
%plot([Obs_max_mid(1) obs_max_loc(1)], [Obs_max_mid(2) obs_max_loc(2)], LineWidth=1, Color='g') % plot h line between obstacle-top and LOS (mid-point)

else

d1      = 0;
d2      = 0;
h       = 0;
F1      = 0;
alph1   = 0;
alph2   = 0;

end

G_info=[d1 d2 h F1 F1_max alph1 alph2]; 

end




