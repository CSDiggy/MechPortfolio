function statedtISS = ISSOrbit(t,stateISS)

xs = stateISS(1);
ys = stateISS(2);
zs = stateISS(3);
xsdot = stateISS(4);
ysdot = stateISS(5);
zsdot = stateISS(6);

% Earth Parameters
Re = 6378e3;
Me = 5.972e24;
G = 6.674e-11;

% Mass of ISS
ISSMass = 450e3;

% Gravitational effects
rISS = stateISS(1:3);
rnormISS = norm(rISS);
rVectorISS = rISS/rnormISS;
gravityForceISS = (-G*Me*ISSMass)/(rnormISS)^2*rVectorISS;

% Kinematics of system
velstateISS = stateISS(4:6);

% Dynamics of system
accelstateISS = gravityForceISS/ISSMass;

statedtISS = [velstateISS;accelstateISS];