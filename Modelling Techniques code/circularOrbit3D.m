function [orbPlot] = circularOrbit3D(heightAboveSurf,incAngle)
%% Plots the projectile motion around the earth in 3D

% This function lets the user plot a 3D circular orbit where the velocity
% varies when inputed with a orbit radius, HEIGHTABOVESURF, and orbit
% inclination angle, INCANGLE.


% Image file for earth render
image_file1 = 'http://upload.wikimedia.org/wikipedia/commons/thumb/c/cd/Land_ocean_ice_2048.jpg/1024px-Land_ocean_ice_2048.jpg';
alpha = 1; % Transparency of render

% Earth parameters
Re = 6378e3;
Rp = 6357e3;
Me = 5.972e24;
G = 6.674e-11;

% Center of ellipsoid coordinates for line tracking
xCenter = 0;
yCenter = 0;
zCenter = 0;

xArray = xCenter(ones(530,1));
yArray = yCenter(ones(530,1));
zArray = zCenter(ones(530,1));

% Initial conditions of projectile
inclineAngle = incAngle;
x0 = 0;
y0 = -Re-heightAboveSurf;
z0 = 0;

% New orbit conditions
xN = 0;
yN = -Re-86e5;
zN = 0;

orbitRadius = norm([x0;y0;z0]);
velCircOrbit = sqrt((G*Me)/orbitRadius);
xVel = velCircOrbit*sind(inclineAngle+22);
yVel = 0;
zVel = -velCircOrbit*cosd(inclineAngle+22);
initialstate1 = [x0;y0;z0;xVel;yVel;zVel];

% Initial conditions of ISS
x0ISS = -Re-400e3;
y0ISS = 0;
z0ISS = 0;

orbitRadiusISS = norm([x0ISS;y0ISS;z0ISS]);
velCircOrbitISS = sqrt((G*Me)/orbitRadiusISS);
xVelISS = 0;
yVelISS = velCircOrbitISS*sind(-30);
zVelISS = -velCircOrbitISS*cosd(-30);

initialstate2 = [x0ISS;y0ISS;z0ISS;xVelISS;yVelISS;zVelISS];

% Initial conditions of elliptical orbit
orbitRadius = norm([x0;y0;z0]);
velCircOrbit = sqrt((G*Me)/orbitRadius);
xVel = velCircOrbit*sind(inclineAngle)*3.2;
yVel = 0;
zVel = -velCircOrbit*cosd(inclineAngle);
initialstate3 = [x0;y0;z0;xVel;yVel;zVel];

% Initial new orbit after transfer
orbitRadiusNew = norm([xN;yN;zN]);
velCircOrbitNew = sqrt((G*Me)/orbitRadiusNew);
xVelNew = velCircOrbitNew*sind(inclineAngle+22);
yVelNew = 0;
zVelNew = -velCircOrbitNew*cosd(inclineAngle+22);
initialstate4 = [xN;yN;zN;xVelNew;yVelNew;zVelNew];

% Time relations
orbitPeriod = sqrt((4*pi^2*orbitRadius^3)/(G*Me));
orbitPeriodISS = sqrt((4*pi^2*orbitRadiusISS^3)/(G*Me));
orbit_numberproj = 4;
orbit_numberISS = 4;
orbit_numberNew = 4;
timespan = [0:30:orbitPeriod*orbit_numberproj*1.2];
timespanISS = [0:15:orbitPeriodISS*orbit_numberISS*1.2];
timespanNewOrb = [0:50:orbitPeriod*orbit_numberNew*1.2];

% Use ODE15s to solve stiff differential equation to high degree of accuracy. Uses
% functions PROJ3DORBIT and ISSORBIT.

% Projectile state 
[tout1,stateout1] = ode15s(@proj3DOrbit,timespan,initialstate1);

% ISS state
[tout2,stateout2] = ode15s(@ISSOrbit,timespanISS,initialstate2);

% Elliptical orbit state
[tout3,stateout3] = ode15s(@proj3DOrbit,timespan,initialstate3);

% New projectile state after transfer
[tout4,stateout4] = ode15s(@proj3DOrbit,timespanNewOrb,initialstate4);

% Recover states in each coordinate direction from STATEOUT1 array for
% projectile.
xstateout = stateout1(:,1);
ystateout = stateout1(:,2);
zstateout = stateout1(:,3);
xVelstateout = stateout1(:,4);
yVelstateout = stateout1(:,5);
zVelstateout = stateout1(:,6);

% Recover states in each coordinate direction from STATEOUT2 array for ISS.
xSstateout = stateout2(:,1);
ySstateout = stateout2(:,2);
zSstateout = stateout2(:,3);

% Recover state in each coordinate direction from STATEOUT3 array for
% elliptical orbit of projectile
xEstateout = stateout3(:,1);
yEstateout = stateout3(:,2);
zEstateout = stateout3(:,3); 

% Recover state in each coordinate direction from STATEOUT4 array for new
% orbit
x2stateout = stateout4(:,1);
y2stateout = stateout4(:,2);
z2stateout = stateout4(:,3);

% Sub plots of displacement
figure(1);
hold on
subplot(5,5,[1,2,3,4,5]);
plot(tout1,xstateout,'r');
xlabel('Time span');
ylabel('X disp');

subplot(5,5,[11,12,13,14,15]);
plot(tout1,ystateout,'r');
xlabel('Time span');
ylabel('Y disp');

subplot(5,5,[21,22,23,24,25]);
plot(tout1,zstateout,'r');
xlabel('Time span');
ylabel('Z disp')
hold off

% Sub plots of velocity
figure(2)
hold on
subplot(5,5,[1,2,3,4,5]);
plot(tout1,xVelstateout,'blue');
xlabel('Time span');
ylabel('X velocity');

subplot(5,5,[11,12,13,14,15]);
plot(tout1,yVelstateout,'blue');
xlabel('Time span');
ylabel('Y velocity');

subplot(5,5,[21,22,23,24,25]);
plot(tout1,zVelstateout,'blue');
xlabel('Time span');
ylabel('Z velocity');
hold off

% Create 'space'
figure('Color','k')

view(0,45);

hold on

set(gca,'Color', 'k', 'Visible','off')
grid on

% Plot box limits set to make rotation origin at centre
xlim([-30e6 30e6])
ylim([-30e6 30e6])
zlim([-30e6 30e6])

% Earth render plot
cdata = imread(image_file1);

[X,Y,Z] = ellipsoid(0, 0, 0, Re, Re, Rp);
orbPlot = surf(X,Y,-Z,'FaceColor','none');
set(orbPlot, 'Facecolor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha, 'Edgecolor', 'none');

% Create fixed global axes
xAx = plot3([0,1.5*Re],[0,0],[0,0],'-.w'); % X axis
xAxT = text(1.5*Re+200,0,0,texlabel('X'),'Color','w');

yAx = plot3([0,0],[0,1.5*Re],[0,0],'-.w'); % Y axis
yAxT = text(0,1.5*Re+200,0,texlabel('Y'),'Color','w');

zAx = plot3([0,0],[0,0],[0,1.5*Re],'-.w'); % Z axis
zAxT = text(0,0,1.5*Re+200,texlabel('Z'),'Color','w');

% Add equitorial line
eqAngles = 1:1:361;
equatorX = (Re*1.0)*cosd(eqAngles);
equatorY = (Re*1.0)*sind(eqAngles);
equatorZ = zeros(1,size(eqAngles,2));

equatorPlot = plot3(equatorX,equatorY,equatorZ,'--w');

plotUpdate = 150;

% Plot orbit lines
ISSorb = plot3(stateout2(:,1),stateout2(:,2),stateout2(:,3),'green','LineWidth',0.5);
projOrbitPlot = plot3(xstateout,ystateout,zstateout,'Color','r','LineWidth',0.5);
ellipProjOrbitplot = plot3(xEstateout,yEstateout,zEstateout,'Color','b','Linewidth',0.5);
newProjOrbitPlot = plot3(x2stateout,y2stateout,z2stateout,'Color','r','LineWidth',1);

% Variable for infinite while loop
whileStat = 0;

% Match Earths tilt
tiltAngle = 23.4;

rotate(orbPlot,[0 1 0],tiltAngle);
rotate(xAx,[0 1 0],tiltAngle);
rotate(xAxT,[0 1 0],tiltAngle);
rotate(yAx,[0 1 0],tiltAngle);
rotate(yAxT,[0 1 0],tiltAngle);
rotate(zAx,[0 1 0],tiltAngle);
rotate(zAxT,[0 1 0],tiltAngle);
rotate(equatorPlot,[0 1 0],tiltAngle);

ISSx = stateout2(:,1);
ISSy = stateout2(:,2);
ISSz = stateout2(:,3);


orbitTimeSpan = 1:1:509;
phProj = plot3(xstateout(1),ystateout(1),zstateout(1),'ro','MarkerSize',8);
phISS = plot3(stateout2(:,1),stateout2(:,2),stateout2(:,3),'Color','green','Marker','o','MarkerSize',8);
phEllipProj = plot3(xEstateout(1),yEstateout(1),zEstateout(1),'bo','MarkerSize',8);
phNewProjOrb = plot3(x2stateout(1),y2stateout(1),z2stateout(1),'ro','MarkerSize',8);

currentRotSpeed=uicontrol(gcf,'Style','text','String','10000','FontSize',12,'BackgroundColor','k','ForegroundColor','w','Units','normalized','Position',[0.69 0.2 0.07 0.05]);
reduceSpeed=uicontrol(gcf,'Style','pushbutton','String','- Speed','Units','normalized','Position',[0.65 0.2 0.05 0.05],'BackgroundColor','k','ForegroundColor','w','Callback',@(~,~) set(currentRotSpeed,'String',num2str(str2num(currentRotSpeed.String)-500)));
increaseSpeed=uicontrol(gcf,'Style','pushbutton','String','+ Speed','Units','normalized','Position',[0.75 0.2 0.05 0.05],'BackgroundColor','k','ForegroundColor','w','Callback',@(~,~) set(currentRotSpeed,'String',num2str(str2num(currentRotSpeed.String)+500)));


% Begin infite while loop
while whileStat == 0
    if str2num(currentRotSpeed.String)>0
        for i=1:length(orbitTimeSpan)

            % Set turnrate dependant to value in text box
            turnRate=360/(86400)/plotUpdate*str2num(currentRotSpeed.String);

            % Update ISS and projectile positions and delete the previous
            % one
            phProj.XData = xstateout(i);
            phProj.YData = ystateout(i);
            phProj.ZData = zstateout(i);

            phISS.XData = ISSx(i);
            phISS.YData = ISSy(i);
            phISS.ZData = ISSz(i);

            phEllipProj.XData = xEstateout(i);
            phEllipProj.YData = yEstateout(i);
            phEllipProj.ZData = zEstateout(i);

            phNewProjOrb.XData = x2stateout(i);
            phNewProjOrb.YData = y2stateout(i);
            phNewProjOrb.ZData = z2stateout(i);
            
            % Connect line from instantaneous positions to centre of earth
            phProjLine = plot3([xArray(i) xstateout(i)],[yArray(i) ystateout(i)],[zArray(i) zstateout(i)],'Color','c','LineWidth',2);
            phISSLine = plot3([xArray(i) xSstateout(i)],[yArray(i) ySstateout(i)],[zArray(i) zSstateout(i)],'Color','c','LineWidth',2);
            %phEllipProj = plot3([xArray(i) xEstateout(i)],[yArray(i) yEstateout(i)],[zArray(i) zEstateout(i)],'Color','c','LineWidth',2);
            phNewProjOrbLine = plot3([xArray(i) x2stateout(i)],[yArray(i) y2stateout(i)],[zArray(i) z2stateout(i)],'Color','c','LineWidth',2);
            
            % Rotate the earth along its origins's Z axis.
            rotate(orbPlot,[0 0 1],turnRate);

            % Rotate earth information with orbPlot
            % Rotate axes and text to stay with globe
            rotate(xAx,[0 0 1],turnRate);
            rotate(xAxT,[0 0 1],turnRate);

            rotate(yAx,[0 0 1],turnRate);
            rotate(yAxT,[0 0 1],turnRate);

            rotate(zAx,[0 0 1],turnRate);
            rotate(zAxT,[0 0 1],turnRate);

            % Rotate equator with globe
            rotate(equatorPlot,[0 0 90],turnRate);

            axis vis3d

            drawnow

            % Delete previous line tracks to keep plot clean
            if i ~= 0
                delete(phISSLine)
                delete(phProjLine)
                delete(phNewProjOrbLine)
            end

            pause(0.00001);
        end
    end
end
