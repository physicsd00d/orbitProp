function statedot = statedot_drag_H(t, state)
% state vector is [rx ry rz vx vy vz]'
mu_earth = 398600.440; % [km^3/sec^2]
% J2 = 1.0826e-3;
% DU = 6378.137;  %km

Cd = 2.3;
mass = 1500;    %[kg]
Area = 20 / (1e3)^2;      %[km^2]

statedot = zeros(6,1);

r = state(1:3);
v = state(4:6);

% Made 1% change to H to evaluate sensitivity
rho0 = 1.225;   %[kg/m^3]
H = 1.01*10e3;     %[m]
DU = 6378.137e3;  %[m]

rho = rho0*exp(-(norm(r)*1000 - DU)/H);  %[kg/m^3]

if (norm(r)*1000 < DU)
    disp('You are underground')
end

%convert rho
rho = rho*(1e3)^3;      %[kg/km^3]

% This is assuming the atmoshpere is at rest wrt inertial space
acc_drag = 0.5*(Cd*Area/mass)*rho*norm(v)*v;

acc = -mu_earth*r/norm(r)^3 - acc_drag;

    statedot(1:3) = state(4:6);
    statedot(4:6) = acc;


end