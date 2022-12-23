%% This script performs a Hohmann Transfer to put the projectile into a higher orbit following its orbital insertion into an altitude of 190km. The manoeuvre transfers the projectile into an orbit of 400km, making it more stable and thus having a much longer decay time.

% Parameters
G = 6.674e-11;
Me = 5.972e24;
projMass = 5.4173e3;

% Important variables
r = [6.568e6, 0, 0];
v = [0, sqrt(G*Me/norm(r(1,:),2)), 0];
t = 0.0; % Start time of simulation
dt = 1; % Time step for data reading
deltaV = 2420; % Change in Velocity required to enter Hohmann Transfer maneouvre
timeToDeltaV = 5431; % Time of applied thrust
mE = 200; % Mass of projectile

% Orbit 1 (Initial)
rplot = r;

while t <= timeToDeltaV

    rMag = norm(r,2);
    accelInstant = -G*Me/rMag^2;
    vNext = v + dt*accelInstant*r/rMag;
    rNext = r+dt*v;

    t = t+dt;

    rplot = [rplot;rNext];
    v = vNext;
    r = rNext;

end

figure(1)
hold on
grid on
title('Trajectory of Hohmann Transfer')
xlim([-5e7 5e7])
ylim([-5e7 5e7])

orbit1 = plot(rplot(:,1), rplot(:,2),'Color','k');

% 2D earth with no details
circAngles = [0:2*pi:500];
earth2DRadius = 6.378e6;
xCirc = earth2DRadius*cosd(circAngles);
yCirc = earth2DRadius*sind(circAngles);
plot(xCirc,yCirc,'Color','blue')

%cometModified(rplot(:,1), rplot(:,2));
axis equal

% Orbit 2 (Transfer)
rplot = r;

Angle = 0;
rPeri = r;

rPeriMag = norm(rPeri,2);
vNext = v + deltaV*[0 1 0];
v = vNext;

while Angle <= 3.1 % A full orbit's angle change is 2pi, half of an orbit is pi.

    rMag = norm(r,2);
    accelInstant = -G*Me/rMag^2;
    vNext = v + dt*accelInstant*r/rMag;
    rNext = r+dt*v;
    t = t+dt;

    rplot = [rplot;rNext];
    v = vNext;
    r = rNext;

    angleDot = dot(rPeri,r);
    cosAngle = angleDot/(rMag*rPeriMag);
    Angle = acos(cosAngle);

end

disp("Event summary:")
disp("Projectile reaches desired orbit height at t = " + t + " seconds.")

orbit2 = plot(rplot(:,1),rplot(:,2),'Color','k');


% Orbit 3 (Final)
rplot = r;
tstart = t;

orbitPeriodHohmann = (2*pi)*sqrt((norm(r,2)).^3/(G*Me));
orbitVelInstant = sqrt((G*Me)/(norm(r,2)));
deltaV = orbitVelInstant - norm(v,2);

vNext = v + deltaV*[0 -1 0];
v = vNext;

while t < (tstart + orbitPeriodHohmann)

    rMag = norm(r,2);
    accelInstant = -G*Me/rMag^2;
    vNext = v + dt*accelInstant*r/rMag;
    rNext = r+dt*v;

    t = t + dt;

    rplot = [rplot;rNext];
    v = vNext;
    r = rNext;

end

orbit3 = plot(rplot(:,1),rplot(:,2),'Color','k');