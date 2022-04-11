clear
format compact

h=0.015;       % timestep of 0.015 time-units
pos=[1 0];     % initial position
vel=[0 0];   % initial velocity

steps=2000;            % number of steps
skip=ceil(steps/100);  % plot only 100 frames

theta=linspace(0,2*pi,100);   % equally spaced points
xc=cos(theta);
yc=sin(theta);

time=0;
for i=1:steps
  x(i)=pos(1);
  y(i)=pos(2);
  time=time+h;
  t(i)=time;
  dist(i)=norm(pos);
  timestep(i)=h;

  E(i)=0.5*norm(vel)^2-1/dist(i)+1;

  % Plot every skip'th frame
  if rem(i,skip)==0
%     plot(xc,yc,'b',x,y,'g-',pos(1),pos(2),'ko',0,0,'ro')
%     title(time/(2*pi))
%     semilogy(dist,timestep,'ro')
    axis equal
    plot(dist,E,'ro');
    drawnow
  end

  % Velocity-Verlet Method
  r=norm(pos);
  accel=-1/r^2 * pos/r;
  pos=pos + h*vel + 0.5*h*h*accel;

  r=norm(pos);
  accelnext=-1/r^2 * pos/r;
  vel=vel + 0.5*h*(accel+accelnext);

  % Variable time-step
  deltav=0.5*h*(accel+accelnext);
  fraction=norm(deltav)/norm(vel);
  if fraction>0.02
    h=h*0.5
  elseif fraction<0.01
    h=h*2
  end

  % Terminate if h too small/large
  if h<1e-6 | h>1000
    break
  end
end

