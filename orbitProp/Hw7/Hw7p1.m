%% Hw7.1
clear all; close all; clc

%% Part d.)

run('/Users/marian/Documents/MATLAB/AA279/Hw7/Gravitation_Parameters.m')

% Equatorial or Mean radii [km]
R_jup = 71492;
R_mar = 3389.9;
R_moo = 1737.53;
R_ven = 6051.8;



disp(['Vmax for Moon     = ' num2str(sqrt(mu_moo/R_moo)) ' km/s'])
disp(['Vmax for Mars     = ' num2str(sqrt(mu_mar/R_mar)) ' km/s'])
disp(['Vmax for Venus    = ' num2str(sqrt(mu_ven/R_ven)) ' km/s'])
disp(['Vmax for Jupiter  = ' num2str(sqrt(mu_jup/R_jup)) ' km/s'])

% Jupiter has the greatest potential delta v change, while the moon has the
% worst.