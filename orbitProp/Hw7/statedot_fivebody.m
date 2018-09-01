function statedot = statedot_fivebody(t, state)
% state vector is [rx ry rz vx vy vz]'
% mu_earth = 1; %DU
% acc_sun = -mu_earth*r/norm(r)^3;
% statedot(1:3) = state(4:6);
% statedot(4:6) = acc;

%Solar System Barycentric Reference Frame (Inertial)

run('/Users/marian/Documents/MATLAB/AA279/Hw7/Gravitation_Parameters.m')

statedot = zeros(30,1);

r_sun = state(1:3);
v_sun = state(4:6);
r_mer = state(7:9);
v_mer = state(10:12);
r_ven = state(13:15);
v_ven = state(16:18);
r_ear = state(19:21);
v_ear = state(22:24);
r_moo = state(25:27);
v_moo = state(28:30);

% account for relative distances!!!

%r_mer_sun is vector from mercury, pointing towards the sun
r_mer_sun = r_sun - r_mer;
r_ven_sun = r_sun - r_ven;
r_ear_sun = r_sun - r_ear;
r_moo_sun = r_sun - r_moo;

acc_sun = -mu_mer*r_mer_sun/norm(r_mer_sun)^3 - mu_ven*r_ven_sun/norm(r_ven_sun)^3 ...
    - mu_ear*r_ear_sun/norm(r_ear_sun)^3 - mu_moo*r_moo_sun/norm(r_moo_sun)^3;


r_sun_mer = r_mer - r_sun;
r_ven_mer = r_mer - r_ven;
r_ear_mer = r_mer - r_ear;
r_moo_mer = r_mer - r_moo;

acc_mer = -mu_sun*r_sun_mer/norm(r_sun_mer)^3 - mu_ven*r_ven_mer/norm(r_ven_mer)^3 ...
    - mu_ear*r_ear_mer/norm(r_ear_mer)^3 - mu_moo*r_moo_mer/norm(r_moo_mer)^3;


r_sun_ven = r_ven - r_sun;
r_mer_ven = r_ven - r_mer;
r_ear_ven = r_ven - r_ear;
r_moo_ven = r_ven - r_moo;

acc_ven = -mu_sun*r_sun_ven/norm(r_sun_ven)^3 - mu_mer*r_mer_ven/norm(r_mer_ven)^3 ...
    - mu_ear*r_ear_ven/norm(r_ear_ven)^3 - mu_moo*r_moo_ven/norm(r_moo_ven)^3;


r_sun_ear = r_ear - r_sun;
r_mer_ear = r_ear - r_mer;
r_ven_ear = r_ear - r_ven;
r_moo_ear = r_ear - r_moo;

acc_ear = -mu_sun*r_sun_ear/norm(r_sun_ear)^3 - mu_mer*r_mer_ear/norm(r_mer_ear)^3 ...
    - mu_ven*r_ven_ear/norm(r_ven_ear)^3 - mu_moo*r_moo_ear/norm(r_moo_ear)^3;


r_sun_moo = r_moo - r_sun;
r_mer_moo = r_moo - r_mer;
r_ven_moo = r_moo - r_ven;
r_ear_moo = r_moo - r_ear;

acc_moo = -mu_sun*r_sun_moo/norm(r_sun_moo)^3 - mu_mer*r_mer_moo/norm(r_mer_moo)^3 ...
    - mu_ven*r_ven_moo/norm(r_ven_moo)^3 - mu_ear*r_ear_moo/norm(r_ear_moo)^3;


statedot(1:3) = v_sun;
statedot(4:6) = acc_sun;
statedot(7:9) = v_mer;
statedot(10:12) = acc_mer;
statedot(13:15) = v_ven;
statedot(16:18) = acc_ven;
statedot(19:21) = v_ear;
statedot(22:24) = acc_ear;
statedot(25:27) = v_moo;
statedot(28:30) = acc_moo;


end