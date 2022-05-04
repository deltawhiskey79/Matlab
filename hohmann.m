%
% Hohmann Transfer Orbit
%

%
% General Parameters
%
%dt = 0.1;
G = 6.67428*10^-11; % Gravity constant
M = 5.9722*10^24; % earth's mass
%T = 20; % duration of circular motion
%t = 0:dt:T; % simulation duration steps
er = 6.371*10^6; % earth's radius 
ar = 2.042*10^6 + er; % parking orbit radius (LEO)
ae = 3.63229*10^8; % perigee distance (earth-moon)
ap = 4.054*10^8; % apogee distance (earth-moon)
a = (ae + ap)/2; % semi major axis of moon's orbit
u = G*M; % u Gravity constant (earth-moon) unit
theta = 0:0.1:360; % one full circle angle @ 0.1 deg

%
% Moon's Orbital Track
%
% Basically any planetary orbit follows elliptical/oval track instead of circular
%
e = 0.0549; % moon's eccentricity
c = e*a;
b = sqrt(a^2 - c^2); % semi-minor axis
lr = b^2/a; % semi-latus rectum or explorer's orbit
r = lr ./ (1 - e*cos(theta*pi/180)); % plotting orbital track
Tm = 2*pi*sqrt(a^3/u); % moon-earth orbit period
dtm = 100; % every time slice is stepped by 100 as the least denominator
tm = Tm*0.98:dtm:Tm; % deliberately set in order to match the orbiter's arrival
mPos = zeros(2, length(tm)); % moon's orbit track data
for i = 1:length(tm)
    mPos(:, i) = (lr ./ (1 - e*cos(2*pi*tm(i)/Tm))) .* [cos(2*pi*tm(i)/Tm); sin(2*pi*tm(i)/Tm)];
end

%
% Parking orbit
%
vo = sqrt(u/ar); % parking orbit velocity
Tp = 2*pi*sqrt(power(ar,3)/u); % LEO period
dtp = 100; % time slice of LEO period stepped by 100 as the least denominator
tp = 0:dtp:Tp; % parking orbit duration
pPos = zeros(2, length(tp)); % parking orbit position data
for i = 1:length(tp)           % Plotting orbiter's track data
    pPos(:,i) = ar * [cos(2*pi*tp(i)/Tp); sin(2*pi*tp(i)/Tp)];
end

%%%%%%%%%%%%%%%%%%%%%%%%
% Goin' to the moon!
%%%%%%%%%%%%%%%%%%%%%%%%

%
% velocity at apogee
%
va = sqrt(u*(2/ar - 1/ap)); % vis viva equation
ao = (ar+ap)/2; % semi-major axis of orbiter
lo = 2*ar*ap/(ap+ar); % semi-latus rectum of orbiter-earth orbit
eo = (ap-ar)/(ap+ar); % orbiter orbit eccentricity
ro = lo ./ (1 - eo*cos(theta*pi/180)); % Orbiter transfer orbit track
To = 2*pi*sqrt(ao^3/u); % Orbiter orbit period
dto = 900;
to = To/2:dto:To; % Assuming the orbiter executes the TLI at perigee
oPos = zeros(2, length(dto)); % Orbiter's track data
%oPos(:,1) = [c-ar; 0];
for i = 1:length(to) % Plotting orbiter's track to intercept the moon
   oPos(:,i) = (lo ./ (1 + eo*cos(2*pi*to(i)/To))) .* [cos(2*pi*to(i)/To); sin(2*pi*to(i)/To)];
end

%
% plot all
%
figure(1);

xo = r .* cos(theta*pi/180)-c; % X axis (elliptic)
yo = r .* sin(theta*pi/180); % Y axis (elliptic)
xe = ar * cos(theta*pi/180)-c; % X axis (parking orbit)
ye = ar * sin(theta*pi/180); % Y axis (parking orbit)
xm = ro .* cos(theta*pi/180)-c; % X axis (orbiter)
ym = ro .* sin(theta*pi/180); % Y axis (orbiter)
plot(xo, yo, '--r'); hold on; % plot track
plot(xe, ye, '--b');
plot(xm, ym, '--bl');
plot(-c, 0 , 'Or');
labels = {'earth'};
title('Hohmann Transfer Orbit (earth-moon)');
text(-c, 0, labels, 'VerticalAlignment','bottom','HorizontalAlignment','right');
  %h1 = plot(pPos(1,1), pPos(1,2), 'xb'); hold on;
h2 = plot(mPos(1,1), mPos(1,2), 'ob');
h3 = plot(oPos(1,1), oPos(1,2), 'xb'); 
% for i = 1:length(oPos)
%    set(h3, 'XData', oPos(1,i)+ap-(ar+c), 'YData', oPos(2,i));
%    pause(0.05);
% end
k = 1; j = 1; % index of loop
while k < length(tm)
%    axis equal;
 if j == length(pPos)
    j = 1; 
 end
 %set(h1, 'XData', pPos(1,j)-c, 'YData', pPos(2,j));
 set(h2, 'XData', mPos(1,k)-c, 'YData', mPos(2,k));
 set(h3, 'XData', oPos(1,k)+ap-(ar+c), 'YData', oPos(2,k));
pause(0.05);
 j = j + 1;
 k = k + 1;
end


%
% Created by: Dimas Widyasastrena, 4th May 2022
%
% "I'm not a genius, no ... seriously; it's just Allah Azza wa Jalla made it easy for
% me"
%