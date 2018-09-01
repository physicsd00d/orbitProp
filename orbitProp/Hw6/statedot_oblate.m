function statedot = statedot_oblate(t, state)
% state vector is [rx ry rz vx vy vz]'
mu_earth = 1; %DU

statedot = zeros(6,1);

r = state(1:3);

acc = -mu_earth*r/norm(r)^3 - 0.2*r/norm(r)^5;

statedot(1:3) = state(4:6);
statedot(4:6) = acc;

end