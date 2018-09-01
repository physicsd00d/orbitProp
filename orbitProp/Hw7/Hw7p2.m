%% Hw7.2
clear all; close all; clc; format long

mu_E = 398600; %km3/s2
R_E = 6378; %km
E_r_M = 384400; %km
mu_M = 4903; %km3/s2
R_M = 1738; %km


%% a.)

v_arrival = [0.7237 0.1872 0]';

h_vec = cross([E_r_M 0 0], v_arrival);

h_mag = norm(h_vec);

disp(['specific angular momentum = ' num2str(h_mag) ' km^2/s'])

%% b.)

r_peri = R_E + 200;

v_peri = h_mag/r_peri;

v_circ = sqrt(mu_E/r_peri);

delta_v_TLI = v_peri - v_circ;

disp(['delta V_TLI = ' num2str(delta_v_TLI) ' km/s'])

%% c.)

% characteristics of transfer ellipse
energy_trans = 0.5*v_peri^2 - mu_E/r_peri;
ecc_trans = sqrt(1+2*energy_trans*h_mag^2/mu_E^2);
p_trans = h_mag^2/mu_E;

% imagine if we had used a hohmann transfer
ecc_hoh = (E_r_M - r_peri)/(E_r_M + r_peri);
p_hoh = 0.5*(E_r_M + r_peri)*(1-ecc_hoh^2);

% enforce constraint r_hoh(nu=180) = r_trans(nu=alpha)
alpha = pi - acos( ((p_trans/p_hoh)*(ecc_hoh-1) + 1)/ecc_trans);

% check that I got the right answer, should = 384400
p_trans/(1 + ecc_trans*cos(alpha));

disp(['alpha = ' num2str(alpha*180/pi) ' deg'])

%% d.)

% alpha_vec = 0:0.001:alpha;
% Y = p_trans./(1+ecc_trans*cos(alpha_vec));
% polar(alpha_vec,Y)
% trapz(alpha_vec,Y)

a_trans = p_trans/(1-ecc_trans^2);

Yfun = @(nu) 0.5*(p_trans./(1+ecc_trans*cos(nu))).^2;
area_swept = quad(Yfun,0,alpha);

area_ellipse = 2*quad(Yfun,0,pi);

tau_ellipse = 2*pi*sqrt(a_trans^3/mu_E);

tau_trans = tau_ellipse*area_swept/area_ellipse;

disp(['transfer time = ' num2str(tau_trans/(3600*24)) ' days'])

%% e.) Just a diagram, see paper

%% f.) 

% Note counterclockwise is positive nu
E_v_s__minM = v_arrival;

% velocity of moon, assuming circular orbit about earth, assuming the
% Earth's position is unaffected by orbit of the moon
E_v_M = [0 sqrt(mu_E/E_r_M) 0]';
% tau_moon = 2*pi*sqrt(E_r_M^3/mu_E)/(3600*24)

M_v_s__mininf = E_v_s__minM - E_v_M;

delta = 2*acos(-M_v_s__mininf(2)/norm(M_v_s__mininf));

ecc_hyperb = 1/sin(delta/2);

disp(['turning angle delta = ' num2str(delta*180/pi) ' degrees'])
disp(['hyperbolic eccentricity = ' num2str(ecc_hyperb) ])

%% g.)

rp_moon = mu_M*(ecc_hyperb-1)/norm(M_v_s__mininf)^2;
b_hyperb = (mu_M/norm(M_v_s__mininf)^2) * sqrt(ecc_hyperb^2-1);

disp(['periselenium = ' num2str(rp_moon) ' km'])
disp(['desired b for lunar insertion = ' num2str(b_hyperb) ' km'])

%% h.)

v_circ_moon = sqrt(mu_M/rp_moon);
v_peri_moon = sqrt( norm(M_v_s__mininf)^2 + 2*mu_M/rp_moon);

delta_v_moon = -v_circ_moon + v_peri_moon;
disp(['delta V for LOI = ' num2str(delta_v_moon) ' km/s'])

%% i.)

M_v_s__plusinf = M_v_s__mininf .* [-1 1 1]';

E_v_s__plusM = M_v_s__plusinf + E_v_M;

% Just using this for graphing purposes
disp(['Velocity of craft on return trip = [' num2str(E_v_s__plusM') '] km/s'])




























