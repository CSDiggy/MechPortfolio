G = 6.674e-11; % Gravitational constant
Me = 5.972e24; % Mass of earth
Re = 6378e3; % Radius of earth
orbitHeight = 190000; % Height of orbit in m
m = mFinal; % Mass of projectile entering orbit
a = (G*Me)/(orbitHeight^2); % Acceleration of gravity at input height, oH
orbitPeriod = sqrt((4*pi^2*orbitHeight^3)/(G*Me)); % Orbital period

% New reading at every point in burn time
dataRead = orbitPeriod/round(burnTime)+1; % Completes circle
motion = -dataRead; % Makes velocity and radius graph cross x = 0. 
for i = 1:1:round(burnTime)

    motion = motion+dataRead; % Ensures that orbit motion starts at t = 0.
    t(i) = motion; 
    orbVelocity(i) = vFinal; 
    rt(i) = orbitHeight; % Creates an array of radii values along with corresponding time values.

end

% Ensure that orbital velocity is constant and stays constant. Graph should be
%straight line
figure
plot(t,orbVelocity);
xlabel('Elapsed time (s)')
ylabel('Velocity (m/s)')

% Ensure that orbital height constant. Graph should be straight line
figure
plot(t,rt);
xlabel('Elapsed time (s)')
ylabel('Orbit Radius (m)')

% Multiple by t to scale 
theta = 2*pi*t/orbitPeriod;
magLine = sqrt((t.^2)+(rt.^2));

figure
polarplot(theta,magLine)

% Could we possibly make decaying orbits?


