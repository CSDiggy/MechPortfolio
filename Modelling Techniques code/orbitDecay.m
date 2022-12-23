function [yDays, yYears] = orbitDecay(aboveSurface, F10, Ap)
% orbitDecayIdea3       Calculates the orbital decay rate from an altitude
% aboveSurface, the solar radio flux, F10 (Between 0 and 300), and geomagnetic index, Ap (0-9). Script made for orbital decay in LEO.
%
% [yDays,yYears] = orbitDecayIdea3(ABOVESURFACE) computes the time taken in
% days for the orbit to fail and commence re-entry. The solution is given
% in days, YDAYS and years, YYEARS.
%
% Function solves for any value of ABOVESURFACE. Main value of ABOVESURFACE
% = 190000 due to assuming circular orbit at apogee in previous scripts.
% Orbit decay takes into consideration the effect of variable molecular
% mass, varying density, solar flux and geomagnetic indexes. Derived
% equations also take into account an approximation of the drag experienced
% by the projectile in LEO in INSTANTDENSITY.

G = 6.674e-11; % Gravitational constant
Me = 5.972e24; % Mass of earth
Re = 6378e3; % Radius of earth
orbitHeight = Re + aboveSurface; 
orbitPeriod = sqrt((4*pi^2*orbitHeight^3)/(G*Me)); % Orbital period
perpA = 20; % Area perpendicular to motion
projMass = 5.3173e3; % Mass while in circular orbit

% Time, T and time increment, INCREMENTT, are in days
incrementT = 0.1;

% Convert to seconds
secondsT = incrementT*3600*24;

% Important decay factors
exoTemp = 900+2.5*(F10-70)+1.5*Ap; % Temperature 
effMolMass = 27-0.012*((aboveSurface/1000)-200); % Effective atmospheric molecular mass. Takes into account variation in molecular mass and compensates for temp variation.
H = exoTemp/effMolMass; % Variable scale height
instantDensity = (6e-10)*(exp(-(((aboveSurface/1000)-175)/H))); % Density at height = aboveSurface

iterCountDec = 0;

while aboveSurface > 180000

    iterCountDec = iterCountDec + incrementT;

    dP = 3*pi*(perpA/projMass)*orbitHeight*instantDensity*secondsT; % Period decay rate due to drag
    orbitHeight = ((orbitPeriod^2*G*Me)/(4*pi^2))^(1/3);
    orbitPeriod = orbitPeriod - dP;

    aboveSurface = orbitHeight - Re;
    
end

reEntryDays = iterCountDec;
reEntryYears = iterCountDec/365;
yYears = round(reEntryYears,4);
yDays = round(reEntryDays,4);

xText = ['The projectile will commence re-entry after ', num2str(yDays), ' days (', num2str(yYears), ' years) if left unattended.'];
disp(xText)







