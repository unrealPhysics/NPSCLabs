clear;
format compact;    % type "help format" to learn what this does

g=9.8;             % acceleration due to gravity
speed=0;          % initial speed is 50 m/s
angle=0;          % beware of radians vs degrees for trig functions

h=0.001;      % timestep of 0.001 second
pos=[0 50];  % start 50 above origin
vel=speed*[cosd(angle) sind(angle)];

rho=1.2; % density of air in kg/m^3
m=45.36; % kg
ironDensity=7800; % kg/m^3
volume=m/ironDensity; % m^3
A=pi*(volume*3/4/pi)^(2/3); % m^2
C_d=0.02;

i=0;
while pos(2)>=0
  i=i+1;
  x(i)=pos(1);
  y(i)=pos(2);

  % Calculating the magnitude of the velocity
  magVel=norm(vel);

  % Calculating the drag force
  dragForce=-0.5*C_d*rho*A*magVel*vel;

  % Calculating the drag acceleration
  dragAccel=dragForce/m;

  % Total acceleration
  accel=dragAccel + [0 -g];

  % Euler method
  pos=pos + h*vel;
  vel=vel + h*accel;
end

% display values

ticsElapsed=i

pos=[0 50];  % start 50 above origin
vel=speed*[cosd(angle) sind(angle)];

m=0.4536; % kg
volume=m/ironDensity; % m^3
A=pi*(volume*3/4/pi)^(2/3); % m^2

i=0;
while i<ticsElapsed
  i=i+1;
  x(i)=pos(1);
  y(i)=pos(2);

  % Calculating the magnitude of the velocity
  magVel=norm(vel);

  % Calculating the drag force
  dragForce=-0.5*C_d*rho*A*magVel*vel;

  % Calculating the drag acceleration
  dragAccel=dragForce/m;

  % Total acceleration
  accel=dragAccel + [0 -g];

  % Euler method
  pos=pos + h*vel;
  vel=vel + h*accel;
end

% display values

disp('Position')
disp(pos)