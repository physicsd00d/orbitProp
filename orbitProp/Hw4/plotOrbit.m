function plotOrbit(orbitalElements, startTime, stopTime, figNum)
%year is assumed to be 2012
%units on elements are assumed to be km, s, rad
%n is assumed to be rev/day
%times are assumed in UTC day + fraction


a = orbitalElements.a;
ecc = orbitalElements.e;
inc = orbitalElements.i;
RAAN = orbitalElements.raan;
omega = orbitalElements.omega;
nu0 = orbitalElements.nu0;
n = orbitalElements.n;
M0 = orbitalElements.M0;
t0 = orbitalElements.t0;

if (ecc == 1)
    error('not designed for parabolic orbits yet')
else
    p = a*(1-ecc^2);
end


% Setup simulation time vector
% t = 1:(24*60*60);
% t = ((to - 105)*24*60*60):(24*60*60);
t = startTime:(1/(24*60*60)):stopTime; 

% % Calculate orbital parameters
% n = (meanmotion*2.0*pi)/(24*60*60);        % [rad/sec]
% p = ((muearth/(n*n))^(1.0/3.0))*(1.0-ecc*ecc); % [km]

eci_C_peri = PERI_C_ECI(RAAN,inc,omega)';             
% Notation note for eci_C_peri: This transforms from perifocal coordinates
% to ECI coordinates according to r_eci = eci_C_peri * r_peri
   
% Set aside space for data
% Nx3 matrices store eci or ecef positions at N datapoints.  Each row
% is one timestep, and columns 1, 2, and 3 hold the three components.
O_r_S__eci  = zeros(length(t),3); % [km]  Sat position in ECI coordinates
O_r_S__ecef = zeros(length(t),3); % [km]  Sat position in ECEF coordinates
lat_S       = zeros(length(t),1); % [rad] Satellite geodetic latitude
lon_S       = zeros(length(t),1); % [rad] Satellite longitude
he_S        = zeros(length(t),1); % [km]  Satellite height above ellipsoid

% Simulate orbit by stepping through time vector
for i = 1:length(t)
%     UTC = UTC_time(1,5,2012,t(i));  %[days]
%      UTC = UTC_time(14,4,2012,t(i));  %[days]
    UTC = t(i);

    % Find mean anomaly and put in range [0, 2*pi)
    MeanAnom = mod(n*(UTC - t0)*2*pi + M0,2*pi);
    
    % Find eccentric anomaly
    EccAnom = EccentricAnomaly(MeanAnom,ecc,1e-10);
     
    % Find position of satellite in perifocal coordinates
    TrueAnom = TrueAnomaly(EccAnom, MeanAnom, ecc);
    r = p./(1+ecc*cos(TrueAnom));
    
    O_r_S__peri = [r*cos(TrueAnom); r*sin(TrueAnom); 0];
    
    % Convert perifocal coordinates to ECI coordinates
    O_r_S__eci(i,:) = eci_C_peri * O_r_S__peri;   % O_r_S__eci(i,:) is a 1x3 row matrix
   
    % Find angle of Greenwich meridian from vernal equinox [0, 2*pi)
%     theta_g = mod(gst2012start + rotEarthDay*(UTC - 1)*2*pi,2*pi);
    theta_g = siderealTime(UTC);
    
    % Convert ECI coordinates to ECEF coordinates
    ecef_C_eci = [cos(theta_g) sin(theta_g) 0; -sin(theta_g) cos(theta_g) 0; 0 0 1];
    
    O_r_S__ecef(i,:) = ecef_C_eci * O_r_S__eci(i,:)';

    
    % Find geodetic latitude, longitude, and height above ellipsoid    
    [lat_S(i), lon_S(i), he_S(i)] = ECEF_To_Geodetic(O_r_S__ecef(i,1), O_r_S__ecef(i,2), O_r_S__ecef(i,3));
end        

% Plot data
% Earth physical constants from Vallado
reearth =   6378.137;               % [km] mean equatorial radius

plotcolor = 'red';
ticki = 60*60; % tick interval - one tickmark every 'ticki' datapoints

% Setup earth plotting data
load('topo.mat', 'topo'); % MATLAB-provided earth data
topoplot = [topo(:,181:360) topo(:,1:180)]; % -180 to +180 longitude
[xeplot, yeplot, zeplot] = ellipsoid(0.0, 0.0, 0.0, ...  % earth sphere
                           reearth, reearth, reearth, 10);

% Plot view from inertial observer
figure(figNum);
clf;
% Draw earth sphere for reference
surface(xeplot, yeplot, zeplot, 'FaceColor', 'blue', 'EdgeColor', 'black');     
hold on;
axis equal;
view(3);
grid on;
set(gca,'XLim', [-50000 +50000], ...
        'YLim', [-50000 +50000], ...
        'ZLim', [-50000 +50000], ...
        'XTick', [-50000:10000:+50000], ...
        'YTick', [-50000:10000:+50000], ...
        'ZTick', [-50000:10000:+50000]);
title('Problem 3.3 - Inertial Observer View');
xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');
plot3(O_r_S__eci(:,1), ...
      O_r_S__eci(:,2), ...
      O_r_S__eci(:,3), ...
      'Color', plotcolor, 'LineStyle', '-', 'LineWidth', 2);
plot3(O_r_S__eci(1:ticki:end,1), ...
      O_r_S__eci(1:ticki:end,2), ...
      O_r_S__eci(1:ticki:end,3), ...
      'Color', plotcolor, 'LineStyle', 'none', ...
      'Marker', 'o', 'MarkerFaceColor', plotcolor, 'MarkerSize', 5);

% Plot view from earth-fixed observer
figure(figNum+1);
clf;
% Draw earth sphere for reference
surface(xeplot, yeplot, zeplot, 'FaceColor', 'blue', 'EdgeColor', 'black');
hold on;
axis equal;
view(3);
grid on;
set(gca,'XLim', [-50000 +50000], ...
        'YLim', [-50000 +50000], ...
        'ZLim', [-50000 +50000], ...
        'XTick', [-50000:10000:+50000], ...
        'YTick', [-50000:10000:+50000], ...
        'ZTick', [-50000:10000:+50000]);
title('Problem 3.3 - Earth-Fixed Observer View');
xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');
plot3(O_r_S__ecef(:,1), ...
      O_r_S__ecef(:,2), ...
      O_r_S__ecef(:,3), ...
      'Color', plotcolor, 'LineStyle', '-', 'LineWidth', 2);
plot3(O_r_S__ecef(1:ticki:end,1), ...
      O_r_S__ecef(1:ticki:end,2), ...
      O_r_S__ecef(1:ticki:end,3), ...
      'Color', plotcolor, 'LineStyle', 'none', ...
      'Marker', 'o', 'MarkerFaceColor', plotcolor, 'MarkerSize', 5);

% Plot groundtrack
figure(figNum+2);
clf;
% Draw earth outline map for reference.  'contour' command may work
% differently in older versions of MATLAB
contour(-180:179, -90:+89, topoplot, [0 0], 'blue');
hold on;
axis equal;
grid on;
set(gca,'XLim', [-180 +180], 'YLim', [-90 +90], ...
        'XTick', [-180:30:+180], 'Ytick', [-90:30:+90]);
title('Problem 3.3 - Ground Track');
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');
% Use 'markers' instead of 'lines' for this plot to avoid distracting
% jumps in plot for data that crosses 180 [deg] longitude
plot(lon_S*180.0/pi, lat_S*180.0/pi, ...
     'Color', plotcolor, 'LineStyle', 'none', ...
     'Marker', 'o', 'MarkerFaceColor', plotcolor, 'MarkerSize', 2); 
plot(lon_S(1:ticki:end)*180.0/pi, lat_S(1:ticki:end)*180.0/pi, ...
     'Color', plotcolor, 'LineStyle', 'none', ...
     'Marker', 'o', 'MarkerFaceColor', plotcolor, 'MarkerSize', 5);



end