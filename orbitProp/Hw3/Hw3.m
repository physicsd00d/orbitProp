%% Problem set 3
clear all; close all; clc;

%% 3.1.a.)

%Stanford is located at
geolat = 37.4226 *pi/180   %[rad]
lon = -122.1654 *pi/180    %[rad]
alt = -9/1E3;   %[km]

[x y z] = Geodetic_To_ECEF(geolat, lon, alt)

%% 3.1.b.)

[gdlat lon alt] = ECEF_To_Geodetic(x, y, z)


