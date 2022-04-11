clear;

global g L         % define g and L as global variables
g=9.8;             % acceleration due to gravity
L=0.2482;          % length in metres (gives period of 1 second)

theta=90*pi/180;   % intial angle (radians)
omega=0;           % intital angular-velocity

x=[theta omega];   % create vector x
h=0.1;            % timestep of 0.1 seconds

% Second-order Runge-Kutta (beta=1)
for i=1:80
  F1=deriv(x);
  F2=deriv(x+h*F1);

  x=x+h/2*(F1+F2);

  % Graphic of the pendulum
  xpend=[0  L*sin(x(1))];
  ypend=[0 -L*cos(x(1))];
  subplot(1,2,1)
  plot(xpend,ypend,'o-')

  axis equal
  axis([-L L -L L])
  thetaVals(i)=x(1);
  time(i)=i*h;
  subplot(1,2,2)
  plot(time,thetaVals,'r-')
  pause(0)
end 
