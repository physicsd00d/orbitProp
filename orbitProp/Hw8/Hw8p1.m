%% Hw8.1
clear all; close all; clc; format long

% Earth physical constants from Vallado
muearth = 398600.4418;              % [km^3/sec^2]
reearth =   6378.137;               % [km] mean equatorial radius
rotEarthRad = 0.0000729211585530;  % [rad/sec]
rotEarthDay = 1.0027379093;  %rev/day
rotEarthDay = rotEarthRad*3600*24/(2*pi);
e2earth =      0.006694385000;      % oblate eccentricity squared

path(path,'Toolbox')
%%

% GPS BIIRM-3 (PRN 12) 
% 1 29601U 06052A   11098.22509229  .00000017  00000-0  10000-3 0  2930 
% 2 29601  55.8420  72.2714 0037332 341.9446  17.9155  2.00551194 32163 
 
yearo_i      = 2011;          % epoch year
to_i         =  98.22509229; % [day] epoch as day of year plus fraction
inc_i        =   55.8420;     % [deg]
RAAN_i       =   72.2714;     % [deg]
ecc_i        =    0.0037332;
omega_i      =  341.9446;     % [deg]
Mo_i         =    17.9155;     % [deg] 
meanmotion_i =    2.00551194; % [rev/day]

%convert the elements from degrees to radians
inc_i = inc_i*pi/180;
RAAN_i = RAAN_i*pi/180;
omega_i = omega_i*pi/180;
Mo_i = Mo_i*pi/180;

n_i = (meanmotion_i*2.0*pi)/(24*60*60);        % [rad/sec]
p_i = ((muearth/(n_i*n_i))^(1.0/3.0))*(1.0-ecc_i*ecc_i); % [km]
a_i = p_i/(1-ecc_i*ecc_i);   %[km] 

tau_i = 2*pi*sqrt(a_i^3/muearth)/(3600*24);

% GPS BIIRM-3 (PRN 12) 
% 1 29601U 06052A   12104.16441319 -.00000023  00000-0  10000-3 0  6951 
% 2 29601  56.1095  57.6374 0039356 355.0140   4.9311  2.00579332 39609 

yearo_f      = 2012;          % epoch year
to_f         =  104.16441319; % [day] epoch as day of year plus fraction
inc_f        =   56.1095;     % [deg]
RAAN_f       =   57.6374;     % [deg]
ecc_f        =    0.0039356;
omega_f      =  355.0140;     % [deg]
Mo_f         =    4.9311;     % [deg] 
meanmotion_f =    2.00579332; % [rev/day]

%convert the elements from degrees to radians
inc_f = inc_f*pi/180;
RAAN_f = RAAN_f*pi/180;
omega_f = omega_f*pi/180;
Mo_f = Mo_f*pi/180;

n_f = (meanmotion_f*2.0*pi)/(24*60*60);        % [rad/sec]
p_f = ((muearth/(n_f*n_f))^(1.0/3.0))*(1.0-ecc_f*ecc_f); % [km]
a_f = p_f/(1-ecc_f*ecc_f);   %[km] 

tau_f = 2*pi*sqrt(a_f^3/muearth)/(3600*24);


end2011 = UTC_time(31,12,2011,24,60,60,0);
start2012 = UTC_time(1,1,2012,0,0,0,0);
delta_time_UTC = (to_f - start2012) + (end2011 - to_i);
delta_time_sec = delta_time_UTC*24*3600;

dRAAN_dt_avg = (RAAN_f - RAAN_i)/delta_time_sec;

J2_est_f = dRAAN_dt_avg * (-2/3) * (p_f/reearth)^2 * 1/(n_f * cos(inc_f))
J2_est_i = dRAAN_dt_avg * (-2/3) * (p_i/reearth)^2 * 1/(n_i * cos(inc_i))

J2_actual = 1.0826e-3






% Not too shabby.  With rounding I get the first two significant digits.













