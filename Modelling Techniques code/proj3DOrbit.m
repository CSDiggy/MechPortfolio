function statedt = proj3DOrbit(t,state)

xs = state(1);
ys = state(2);
zs = state(3);
xsdot = state(4);
ysdot = state(5);
zsdot = state(6);

% Earth parameter
Re = 6378e3;
Me = 5.972e24;
G = 6.674e-11;

% Projectile mass in orbit
projMass = 5.3173e3;

% Gravitational effects
r = state(1:3);
rnorm = norm(r);
rVector = r/rnorm;
gravityForce = (-G*Me*projMass)/(rnorm)^2*rVector;

% Kinematics of system
velstate = state(4:6);

% Dynamics of system
accelstate = gravityForce/projMass;

statedt = [velstate;accelstate];