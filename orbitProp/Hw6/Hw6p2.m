%% Hw6
format long

%% 6.2.a)
clear all; close all; clc;

tau_gs = 86164; %s      geostationary period
mu_earth = 398600.4418; %km^3/s^2
DU = 6378.137;

a_gs = ((tau_gs/(2*pi))^2 * mu_earth)^(1/3);

rA = 200 + DU;
rB = a_gs;

a_t = 0.5*(rA + rB);

vA_minus = sqrt(mu_earth/rA);
vB_plus = sqrt(mu_earth/rB);

vA_plus = sqrt(2*mu_earth*(1/rA - 1/(rA+rB)));
vB_minus = sqrt(2*mu_earth*(1/rB - 1/(rA+rB)));

delta_vA = vA_plus - vA_minus;
delta_vB = vB_plus - vB_minus;
delta_v = delta_vA + delta_vB
tau_transfer = pi*sqrt(a_t^3/mu_earth)

disp(['delta_v = ' num2str(delta_v) ' km/s'])
disp(['transfer time = ' num2str(tau_transfer/3600) ' hours']) 


%% 6.2.b)
clear all;

tau_gs = 86164; %s      geostationary period
mu_earth = 398600.4418; %km^3/s^2
DU = 6378.137;

rA = 200 + DU;
rB = 384400;

a_t = 0.5*(rA + rB);

vA_minus = sqrt(mu_earth/rA);
vA_plus = sqrt(2*mu_earth*(1/rA - 1/(rA+rB)));
delta_vA = vA_plus - vA_minus;
disp(['delta_vA = ' num2str(delta_vA) ' km/s'])

%% 6.2.c)

% No, I am not surprised because the calculation to get to the moon is just
% a collision course.  I would be more surprised if the delta-V to park in
% an orbit at the moon were less than for going to geostationary because
% those are more equivalent scenarios to compare.

% Also, there isn't really much in space to oppose the motion of a moving
% body, so a little bit of thrust can go a long way.















