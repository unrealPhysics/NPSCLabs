clear;

h=0.05;       % timestep of 0.05 time-units
pos=[1 0];    % initial position
vel=[0 1.1];  % initial velocity

steps=100;              % number of steps
skip=ceil(steps/100);   % plot only 100 frames

theta=linspace(0,2*pi,100);   % equally spaced points
xc=cos(theta);
yc=sin(theta);

for i=1:steps
  x(i)=pos(1);
  y(i)=pos(2);
  time=(i-1)*h/(2*pi);
  t(i)=time;
  magVel(i)=norm(vel);
  magPos(i)=norm(pos);

  if rem(i,skip)==0
      subplot(2,1,1);
      plot(t,magVel,'b')

      subplot(2,1,2);
      plot(t,magPos,'b')
  end

  r=norm(pos);
  accel=-1/r^2 * pos/r;

  vel=vel+h*accel;
 
  % Verlet Method
  if i==1
    next=pos + h*vel + 0.5*h*h*accel;
  else
    next=2*pos - prev + h*h*accel;
  end
 
  prev=pos;
  pos=next;
end

disp(max(magVel)/min(magVel))
disp(max(magPos)/min(magPos))

disp(pos)
