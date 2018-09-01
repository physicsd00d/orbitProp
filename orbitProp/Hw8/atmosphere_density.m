function rho = atmosphere_density(r)
% r is in [m], answer is in kg/m^3

rho0 = 1.225;   %[kg/m^3]
H = 10e3;     %[m]
DU = 6378.137e3;  %[m]

rho = rho0*exp(-(r - DU)/H);

if (r < DU)
    disp('You are underground')
end

end