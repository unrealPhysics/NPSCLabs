clear;

h=0.000001;       % timestep of 0.05 time-units
pos=[1 0 0];    % initial position
vel=[0 1.1 0];  % initial velocity

steps=25000000;              % number of steps
skip=ceil(steps/100);   % plot only 100 frames

theta=linspace(0,2*pi,100);   % equally spaced points
xc=cos(theta);
yc=sin(theta);

for i=1:steps
  x(i)=pos(1);
  y(i)=pos(2);
  time=(i-1)*h/(2*pi);

  L=cross(pos,vel);
  angularMomentum(i)=L(3);
  t(i)=time;

%   if rem(i,skip)==0
%     plot(xc,yc,'b',x,y,'g-',pos(1),pos(2),'ko',0,0,'ro')
%     title(time)
%     axis equal;
%     pause(0);
%   end

  if rem(i,skip)==0
    plot(t,angularMomentum-1.1,'g-')
    pause(0);
  end

  r=norm(pos);
  accel=-1/r^2 * pos/r;
 
  % Velocity Verlet Method
  pos=pos + h*vel + 0.5*h^2*accel;
  r=norm(pos);
  accelNext=-1/r^2 * pos/r;
  vel=vel + 0.5*h*(accel+accelNext);

end
