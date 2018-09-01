function [gdlat lon alt] = ECEF_To_Geodetic(x, y, z)


R_equator = 6378.137;   %[km]
ecc_Earth = 0.081819221456;


rdelta = sqrt(x^2 + y^2);   

tol = 1e-10;

gdlat = asin(z/sqrt(x^2 + y^2 + z^2));
err = 50;

if ((abs(x) <= tol) && (abs(y) <= tol))
    lon = 0;
else
    lon = atan2(y,x);    %[deg]
end


while( err > tol )
    Nlat = R_equator/sqrt(1-(ecc_Earth*sin(gdlat))^2);

    gdlat_new = atan((z+Nlat*sin(gdlat)*ecc_Earth^2)/rdelta);

    err = abs(gdlat_new - gdlat);
    gdlat = gdlat_new;
end

alt = rdelta/cos(gdlat) - Nlat;



end