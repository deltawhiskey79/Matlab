%
% General parameters
%
r = 50;
dt = 0.1;
T = 20; % period of circular motion
t = 0:dt:T; % simulation duration
mPos = zeros(2, size(t,2));
mPos(:,1) = [r; 0]; % initial position of circular orbit object 
ePos = zeros(2, size(t,2));
ePos(:,1) = [a; 0]; % initial position of elliptic orbit object

%
% Circular orbit track
%
theta = 0:0.1:360;
x = r * cos(theta*pi/180);
y = r * sin(theta*pi/180);

%
% Elliptic orbit track
%
a = 8; % semimajor axis
b = 3; % semiminor axis
c = sqrt(a^2 - b^2); % focal
e = sqrt(1 - b^2/a^2); % eccentricity
l = a*(1 - e^2); % latus rectum
r1 = l ./ (1 + e*cos(theta*pi/180));
x1 = r1 .* cos(theta*pi/180) + c; % puts the ellipse at the center
y1 = r1 .* sin(theta*pi/180);

% Simulation data of circual orbit 
for i = 2:length(t)   
    mPos(:, i) = r*[cos(2*pi*t(i)/T); sin(2*pi*t(i)/T)];
end

% Simulation data of elliptical orbit
for i = 2:length(t)
    ePos(:, i) = (l./(1 + e*cos(2*pi*t(i)/T))).*[cos(2*pi*t(i)/T); sin(2*pi*t(i)/T)];
end

%plot circular orbit simulation
axis equal;
figure(1)
h1 = plot(mPos(1, 1), mPos(1, 2), 'ob');hold on;
plot(x, y, '--r');
plot(0, 0, 'Or');
for i = 1:length(mPos)
   axis equal;
   set(h1, 'XData', mPos(1,i), 'YData', mPos(2,i)); 
   pause(0.05);
   xlim([-60 60]);
   ylim([-60 60]);
end

%plot ellipse orbit simulation
axis equal;
figure(2)
h2 = plot(ePos(1, 1), ePos(1, 2), 'ob');hold on;
plot(x1,y1,'--r');
plot(-c,0, 'Ob');
for i = 1:length(ePos)
    axis equal;
    set(h2, 'XData', ePos(1,i)+c, 'Ydata', ePos(2,i));
    pause(0.05);
    xlim([-10 10]);
    ylim([-6 6]);
end

%
% Created by : Dimas Widyasastrena @2022
%

