%
% Hohmann Transfer Orbit
%

%
% General Parameters
%
dt = 0.1;
T = 20; % duration of circular motion
t = 0:dt:T; % simulation duration steps
ar = 3*10^6; % parking orbit radius (LEO) in km
u = 4.03475*10^14; % u Gravity constant (earth-moon)
theta = 0:0.1:360; % one full circle angle @ 0.1 deg

%
% Orbital Track
%
% Basically all orbits follow elliptical/oval track instead of circular
%
e = 0.0549; % moon's eccentricity
ap = 3.844*10^5; % apogee point from earth (explorer intercept)
c = e*a;
b = sqrt(a^2 - c^2);
lr = b^2/a; % semi-latus rectum or explorer's orbit
a = (ap + ar)/2; % semi major axis
r = lr ./ (1 - e*cos(theta*pi/180));

%
% Parking orbit
%
vo = sqrt(u/ar); % parking orbit velocity

%
% velocity at apogee
%
va = sqrt(u*(2/ar - 1/a)); % vis viva equation




