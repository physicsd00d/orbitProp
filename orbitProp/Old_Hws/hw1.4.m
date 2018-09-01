%% HW1 Prob 4
clear all; close all; clc;

%% a.)
disp('part a')
clear all; close all;
r = [0; 0.9; 0.3];
v = [1.1; 0; 0.2];
en = norm(v)^2/2 - 1/norm(r)

th = acosd(r'*v/(norm(r)*norm(v)))
hmag = norm(r)*norm(v)*sind(th)

ecc = sqrt(1+2*en*hmag^2)

% thus ellipse

%% b.)
disp('part b')
clear all; close all;

rmag = 2;
vmag = 1.5;

en = vmag^2/2 - 1/rmag

% thus hyperbola


%% c.)
disp('part c')
clear all; close all;

rp_mag = 2;
hmag = 2;

vp_mag = hmag/rp_mag

en = vp_mag^2/2 - 1/rp_mag

% thus parabola

%% d.)
disp('part d')
clear all; close all;

r = [0; 0; 0.4];
v = [1.5; 0; 1.0];
en = norm(v)^2/2 - 1/norm(r)

th = acosd(r'*v/(norm(r)*norm(v)))
hmag = norm(r)*norm(v)*sind(th)

ecc = sqrt(1+2*en*hmag^2)

% thus ellipse again

%% e.)
disp('part e')
clear all; close all;

en = -1/7;
p = 3.5;

hmag = sqrt(p)
ecc = sqrt(1+2*en*hmag^2)










