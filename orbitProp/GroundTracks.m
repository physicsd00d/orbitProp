clear all; close all; clc;

DU = 6378.1;    %km 
run('/Users/marian/Documents/MATLAB/plot_trajectory.m')
%% Plot the nominal trajectory on the Earth
% http://www.mathworks.com/products/matlab/demos.html?file=/products/demos/shipping/matlab/earthmap.html

% load('topo.mat','topo','topomap1');
% whos topo topomap1
% 
% colormap(topomap1);
% 
% % Create the surface.
% [xe,ye,ze] = sphere(50);
% props.AmbientStrength = 0.1;
% props.DiffuseStrength = 1;
% props.SpecularColorReflectance = .5;
% props.SpecularExponent = 20;
% props.SpecularStrength = 1;
% props.FaceColor= 'texture';
% props.EdgeColor = 'none';
% props.FaceLighting = 'phong';
% props.Cdata = topo;
% surface(-xe*DU,-ye*DU,ze*DU,props); %negative signs so that their coord convention aligns with mine
% 
% 
% % Set the view.
% axis square off
% axis equal
% view(3)


%%

load('topo.mat','topo','topomap1');
whos topo topomap1


scrsz = get(0,'ScreenSize');
fighandle = figure('Position',[1 scrsz(4)/2 scrsz(3)/1.1 scrsz(4)/1.5])

reearth = DU;
% Setup earth plotting data
% load('topo.mat', 'topo'); % MATLAB-provided earth data
topoplot = [topo(:,181:360) topo(:,1:180)]; % -180 to +180 longitude
% [xeplot, yeplot, zeplot] = ellipsoid(0.0, 0.0, 0.0, ...  % earth sphere
%                            reearth, reearth, reearth, 10);

contour(-180:179, -90:+89, topoplot, [0 0], 'blue');

hold on


image([0 180],[-90 90],topo(:,1:180),'CDataMapping', 'scaled');
hold on 
image([-180 0],[-90 90],topo(:,181:360),'CDataMapping', 'scaled');
colormap(topomap1);

axis equal
box on
set(gca,'XLim',[-180 180],'YLim',[-90 90], ...
    'XTick',[-180 -120 -60 0 60 120 180], ...
    'Ytick',[-90 -60 -30 0 30 60 90]);

for i=1:length(O_r_S__geodetic)
    hold on
    plot(O_r_S__geodetic{i}(:,2)*180/pi, O_r_S__geodetic{i}(:,1)*180/pi,'r','LineWidth',1)
end




















