%% Hw6.4
clear all; close all; clc;

path(path,'Toolbox6')
%% Part b.)
r0 = [4 0 0]';  %DU
v0 = [0 0.5 0]';  %DU/TU
state0 = [r0; v0];

[a ecc inc raan aop nu0 meanmotion M0] = getOrbitalElements(r0, v0);

tau = 2*pi*sqrt(a^3);   %orbit period

t_vec = 0:0.5:tau;


options = odeset('RelTol', 1e-3, 'AbsTol', 1e-6);
[t_out, state] = ode45(@statedot, t_vec, state0, options) ;

figure(1)
plot(state(:,1),state(:,2))
axis([-6 6 -6 6])
axis equal
title('Circular Orbit over single period with ode45')
xlabel('x [DU]')
ylabel('y [DU]')


%% Part c.)

t_vec = 0:0.5:10*tau;

options = odeset('RelTol', 1e-3, 'AbsTol', 1e-6);
[t_out, state] = ode45(@statedot, t_vec, state0, options) ;

figure(2)
plot(state(:,1),state(:,2))
axis([-6 6 -6 6])
axis equal
title('Circular Orbit over 10 periods with ode45')
xlabel('x [DU]')
ylabel('y [DU]')


% The slow change in the orbit's shape is definitely a numerical effect.
% Nowhere have I accounted for any forces other than simple center-pointing
% gravity.  I find it especially weird that when you nail the period, the
% orbit hooks up with itself and makes a closed circle, but when you
% over-estimate the period, the orbit no longer closes back upon itself.  I
% would have thought that the error there would be roughly the same, but I
% guess the ode45 takes future guesses into account so if there is no
% future, you get a different answer?


%% Part d.)

t_vec = 0:0.5:10*tau;

options = odeset('RelTol', 1e-3, 'AbsTol', 1e-6);
[t_out, state] = ode113(@statedot, t_vec, state0, options) ;

figure(3)
plot(state(:,1),state(:,2))
axis([-6 6 -6 6])
axis equal
title('Circular Orbit over 10 periods with ode113')
xlabel('x [DU]')
ylabel('y [DU]')


% This method is doing WAY better than ode45.  You can still see some
% numerical error, but it's much smaller.


%% Part e.)

r0 = [4 0 0]';  %DU
v0 = [0 0.4 0]';  %DU/TU
state0 = [r0; v0];

[a ecc inc raan aop nu0 meanmotion M0] = getOrbitalElements(r0, v0);

tau = 2*pi*sqrt(a^3);   %orbit period

t_vec = 0:0.5:10*tau;

options = odeset('RelTol', 1e-3, 'AbsTol', 1e-6);
[t_out, state] = ode113(@statedot, t_vec, state0, options) ;

figure(4)
plot(state(:,1),state(:,2))
axis([-6 6 -6 6])
axis equal
title('Elliptical Orbit over 10 periods with ode113')
xlabel('x [DU]')
ylabel('y [DU]')


%% Part f.)

options = odeset('RelTol', 1e-3, 'AbsTol', 1e-6);
[t_out, state] = ode113(@statedot_oblate, t_vec, state0, options) ;

figure(5)
plot(state(:,1),state(:,2))
axis([-6 6 -6 6])
axis equal
title('Elliptical Orbit around Oblate Earth over 10 periods with ode113')
xlabel('x [DU]')
ylabel('y [DU]')

% The orbit changes here are still SOMEWHAT due to numerical error, but
% they are very small compared to the dominant effect of Earth's
% oblateness.  What we are seeing here is a real effect.









