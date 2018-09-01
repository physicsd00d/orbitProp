function [ EccAnom ] = EccentricAnomaly( MeanAnom, ecc, tol )

if (tol <= 0)
    error('Your tolerance is a little unrealistic, huh?')
end

len = length(MeanAnom);
EccAnom = zeros(len,1);

for i=1:len
    En = 0;
    Mn = 0;
    err = 1E5;
    
    while (err >= tol)
        EccAnom(i) = En + (MeanAnom(i) - En + ecc*sin(En))/(1 - ecc*cos(En));
        err = abs(EccAnom(i)-En);
        En = EccAnom(i);
    end
end
