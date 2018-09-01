function [x y z] = Geodetic_To_ECEF(gdlat, lon, alt)
% Wants the angles in degrees

R_equator = 6378.137;   %[km]
ecc_Earth = 0.081819221456;

Nlat = R_equator/sqrt(1-(ecc_Earth*sin(gdlat))^2);

x = (Nlat + alt)*cos(gdlat)*cos(lon);
y = (Nlat + alt)*cos(gdlat)*sin(lon);
z = (Nlat*(1-ecc_Earth^2) + alt)*sin(gdlat);
end