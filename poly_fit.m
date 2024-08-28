%==========================================================================
% Dr. Muneer Al-Zubi
% Linkdin: https://www.linkedin.com/in/muneeralzubi85/
% Email:   muneermaz85@gmail.com

% For more details, the reader can refer to the following article:
% M M. Al-Zubi; M.-S. Alouini, "End-to-End Modelling and Simulation of NLOS Sub-6 GHz Backhaul via Diffraction for
% Internet Connectivity of Rural Areas", IEEE Open Journal of the Communications Society, 2023.
% https://ieeexplore.ieee.org/document/10306284
%==========================================================================

function [x0, y0, x_a, z_a, x_b, z_b] = poly_fit (x, z, tx_loc,rx_loc)

% x is points describe the peak interval of terrain and y is its elevation values.
% p =[p2 p1 p0] is parabola coefficent where z = p1 * x^2 + p2 * x + p3 
% tx_loc/rx_loc is location above seal level

p        =   polyfit(x,z,2); % coefficents of poly function (parabola) where  parabola equation : z = a(x-x0)^2 + y0 = (a) * x^2 - (2*a*x0) * x + (y0+a*x0^2)
                             % Vertex of parabola = (x0, y0)
a        =   p(1); 
x0       =   - p(2) / (2*a); 
y0       =   p(3) - a * x0^2;

%===============================
% tx 
syms xi
m1          =  2 * p(1) * xi + p(2) ; % slope of tangent line equal derivatve of  z = p1 * x^2 + p2 * x + p3
m2          = ( ( p(1) * xi.^2 + p(2) * xi + p(3)) - tx_loc(2) ) ./(xi - tx_loc(1) ) ; % slop equal to diff(y)/diff(x)

x_a = solve(m1-m2==0);
z_a  = p(1) * x_a.^2 + p(2) * x_a + p(3);

 if (z_a(1)>=0)
 z_a = double(x_a(1));
 x_a = double(x_a(1));
 else
 z_a = double(z_a(2));
 x_a = double(x_a(2));
 end
 
%===============================
% rx
syms xi
m1          =  2 * p(1) * xi + p(2) ; % slope of tangent line equal derivatve of  z = p1 * x^2 + p2 * x + p3
m2          = ( ( p(1) * xi.^2 + p(2) * xi + p(3)) - rx_loc(2) ) ./(xi - rx_loc(1) ) ; % slop equal to diff(y)/diff(x)

x_b  = solve(m1-m2==0);
z_b  = p(1) * x_b.^2 + p(2) * x_b + p(3);

 if (z_b(1)>=0)
 z_b = double(z_b(1));
 x_b = double(x_b(1));
 else
 z_b = double(z_b(2));
 x_b = double(x_b(2));
 end
%===============================

z1       =   polyval(p,x); % get y values of obtained parabola (z1 = a*(x-x0).^2 + y0;)

hold on;
plot(x, z1,'r', LineWidth=2)
plot(x0,y0,'xg', LineWidth=2)
plot([tx_loc(1) x_a], [tx_loc(2)  z_a],'r', LineWidth=2)
plot([rx_loc(1) x_b], [rx_loc(2)  z_b],'r', LineWidth=2)

end