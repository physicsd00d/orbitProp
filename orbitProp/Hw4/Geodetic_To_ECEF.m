function ECEF = Geodetic_To_ECEF(gdlat, lon, alt)

R_equator = 6378.137;   %[km]
ecc_Earth = 0.081819221456;

Nlat = R_equator/sqrt(1-(ecc_Earth*sin(gdlat))^2);

ECEF = [(Nlat + alt)*cos(gdlat)*cos(lon);
        (Nlat + alt)*cos(gdlat)*sin(lon);
        (Nlat*(1-ecc_Earth^2) + alt)*sin(gdlat)];


end