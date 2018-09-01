function RotMat = PERI_C_ECI(Rasc, Inc, AoP)
%radians!

RotMat = ...
    [cos(Rasc)*cos(AoP) - cos(Inc)*sin(Rasc)*sin(AoP), cos(AoP)*sin(Rasc) + cos(Rasc)*cos(Inc)*sin(AoP), sin(Inc)*sin(AoP);...
    -cos(Rasc)*sin(AoP) - cos(Inc)*cos(AoP)*sin(Rasc), cos(Rasc)*cos(Inc)*cos(AoP) - sin(Rasc)*sin(AoP), cos(AoP)*sin(Inc);...
     sin(Rasc)*sin(Inc)                                 , -cos(Rasc)*sin(Inc)                                , cos(Inc)];




end