function TrueAnom = TrueAnomaly(EccAnom, MeanAnom, ecc)

num = cos(EccAnom) - ecc;
den = 1 - ecc*cos(EccAnom);

half_plane = floor(MeanAnom/pi);  %which half plane are we in?

TrueAnom = half_plane*pi + acos((-1).^mod(half_plane,2) .* num./den);

TrueAnom = mod(TrueAnom,2*pi);




end