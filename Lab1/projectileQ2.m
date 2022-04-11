clear;
format compact;    % type "help format" to learn what this does

g=9.8;             % acceleration due to gravity
speed=50;          % initial speed is 50 m/s
angle=30;          % beware of radians vs degrees for trig functions

h=0.001;      % timestep of 0.001 second
pos=[0 0];  % start at the origin
vel=speed*[cosd(angle) sind(angle)];

rho=1.2; % density of air in kg/m^3
m=0.145; % kg
r=0.037; % m
A=pi*r^2; % m^2
C_d=0.35;

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

% plot the trajectory as points connected by lines
plot(x,y,'o-')
xlabel('distance (m)')
ylabel('height (m)')

% linear interpolation to estimate the range of the projectile
slope=(pos(2)-y(end))/(pos(1)-x(end));
range=x(end)-y(end)/slope;
disp('Range')
disp(range)
%disp('True Range')
%disp(speed^2*sind(2*angle)/g)