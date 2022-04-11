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
  F2=deriv(x+0.5*h*F1);
  F3=deriv(x+0.5*h*F2);
  F4=deriv(x+h*F3);

  x=x+h/6*(F1+2*F2+2*F3+F4);

%   % Velocity-Verlet Method
%   alpha=-g/L*sin(x(1));
%   x(1)=x(1) + h*x(2) + 0.5*h*h*alpha;
%   alphaNext=-g/L*sin(x(1));
%   x(2)=x(2) + 0.5*h*(alpha+alphaNext);

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
