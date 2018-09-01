function statedot = statedot_oblate(t, state)
% state vector is [rx ry rz vx vy vz]'
mu_earth = 398600.440; % [km^3/sec^2]
J2 = 1.0826e-3;
DU = 6378.137;  %km


statedot = zeros(6,1);

r = state(1:3);
rx = r(1);
ry = r(2);
rz = r(3);

k_hat = [0 0 1]';

acc_oblate = (mu_earth*J2*(DU^2)/2) * ...
    ( (6*rz/norm(r)^5)*k_hat + (3/norm(r)^5 - 15*(rz^2/norm(r)^7))*r );

acc = -mu_earth*r/norm(r)^3 - acc_oblate;

statedot(1:3) = state(4:6);
statedot(4:6) = acc;

end