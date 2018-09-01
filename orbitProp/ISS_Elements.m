% Orbital Elements of ISS for Francisco at time of STS-111 launch

path(path,'Toolbox')

AU = 149597870.691;  %km
DU = 6378.137;  %km
TU = 806.80415; %sec
VU = DU/TU;

% at this instant, sidereal time is about 14.18hrs
thetag = 14.18*2*pi/24;

r = [2.844160626273773E-05 -3.384067431755940E-05 -9.565249794839739E-06]'*AU;
v = [1.083289994539589E-03  1.982654671773565E-03 -3.813957210308104E-03]'*AU/(24*3600);

% go canonical
[a ecc inc raan aop nu0 meanmotion M0] = getOrbitalElements(r/DU,v/VU)

% NOTE: the coordinate system for r and v is weird...not what it should be
% at all.  the elements output here are definitely wrong :(