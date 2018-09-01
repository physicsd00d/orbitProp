function theta_g = siderealTime(UTC)

% Earth physical constants from Vallado
rotEarthRad = 0.0000729211585530;  % [rad/sec]
rotEarthDay = 1.0027379093;  %rev/day
rotEarthDay = rotEarthRad*3600*24/(2*pi);

% Greenwich Sidereal Time at 0000h 1 January 2012 UTC
% from US Naval Observatory website
gst2012start = 6.6706801; % [sidereal hours] from vernal equinox

%convert to rad
gst2012start = gst2012start * 2*pi/24;    %[rad]

% Find angle of Greenwich meridian from vernal equinox [0, 2*pi)
% theta_g = mod(gst2012start + rotEarthDay*(UTC - 1)*2*pi,2*pi);
theta_g = mod(gst2012start + rotEarthRad*(UTC - 1)*3600*24,2*pi);

end