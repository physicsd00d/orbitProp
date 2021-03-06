%% HW1 Prob 4
clear all; close all; clc;

%% a.)
p = 1;
ecc = 0.7;

nu = 1:360;
r = p./(1 + ecc*cosd(nu));

rp = r(end);
ra = r(180);

a = 0.5*(rp+ra);

x = r.*cosd(nu);
y = r.*sind(nu);

figure(1)
plot(x,y)
hold on
plot([0 0],[0 p],'black')
hold on
plot([x(end) (x(end)-a)], [0 0], 'g')
hold on
plot([0 x(end)],[0 y(end)],'r','LineWidth',2)
axis('equal')
title(['a.) ecc=' num2str(ecc) '  p=1  a=' num2str(a) '  rp=' num2str(rp)])

figure(5)
set(gcf, 'Visible', 'off') 
plot(x,y,'-','LineWidth',2)
legend('0.7')


%% b.)
clear all;

p = 1;
ecc = 0.0;

nu = 1:360;
r = p./(1 + ecc*cosd(nu));

rp = r(end);
ra = r(180);

a = 0.5*(rp+ra);

x = r.*cosd(nu);
y = r.*sind(nu);

figure(2)
plot(x,y)
hold on
plot([0 0],[0 p],'black')
hold on
plot([x(end) (x(end)-a)], [0 0], 'g')
hold on
plot([0 x(end)],[0 y(end)],'r','LineWidth',2)
axis('equal')
title(['b.) ecc=' num2str(ecc) '  p=1  a=' num2str(a) '  rp=' num2str(rp)])

figure(5)
set(gcf, 'Visible', 'off') 
hold on
plot(x,y,'-.','LineWidth',2)
legend('0.7','0.0')

%% c.)
clear all;

p = 1;
ecc = 2.2;

%turning angle
delta = 2*asind(1/ecc);
nu_limit = (95);   %actually equal to 90 + delta/2, but that's too big to plot nicely

th = (90 + delta/2)*ones(10,1);

m = tand(90+delta/2);






nu = [1:nu_limit NaN (360-nu_limit):360];
r = p./(1 + ecc*cosd(nu));

rp = p/(1+ecc);
% ra = r(180);

a = rp/(1-ecc);


x_asym = -0.2:0.05:(rp-a);
y_asym = m*x_asym - m*(rp-a);




x = r.*cosd(nu);
y = r.*sind(nu);

figure(3)
plot(x,y)
hold on
plot([0 0],[0 p],'black')   %p
hold on
plot([x(end) (x(end)-a)], [0 0], 'g')   %a
hold on
plot([0 x(end)],[0 y(end)],'r','LineWidth',2)   %rp
axis('equal')
hold on
plot(x_asym,y_asym,'--')  %asymptote above
hold on
plot(x_asym,-y_asym,'--') %asymptote below
title(['c.) ecc=' num2str(ecc) '  p=1  a=' num2str(a) '  rp=' num2str(rp)])

figure(5)
set(gcf, 'Visible', 'off') 
hold on
plot(x,y,'--','LineWidth',2)

%% d.)
clear all;

p = 1;
ecc = 1.0;

nu_limit = 110;
nu = [1:nu_limit NaN (360-nu_limit):360];
r = p./(1 + ecc*cosd(nu));

rp = r(end);
% ra = r(180);

a = NaN;

x = r.*cosd(nu);
y = r.*sind(nu);

figure(4)
plot(x,y)
hold on
plot([0 0],[0 p],'black')
hold on
plot([x(end) (x(end)-a)], [0 0], 'g')
hold on
plot([0 x(end)],[0 y(end)],'r','LineWidth',2)
axis('equal')
title(['d.) ecc=' num2str(ecc) '  p=1  a=' num2str(a) '  rp=' num2str(rp)]) 

figure(5)
set(gcf, 'Visible', 'off') 
hold on
plot(x,y,':','LineWidth',2)


%% e.)

figure(5)
hold on
scatter(0,0,'black','filled')
axis('equal')
legend('0.7','0.0','2.2','1.0')
title('e.) Composite Plot')

























