clear;

g=9.8;             % acceleration due to gravity
L=0.2482;          % length in metres (gives period of 1 second)

theta=10*pi/180;   % initial angle (in radians)
omega=0;           % initial angular-velocity (in radians/second)
h=0.01;            % timestep (in seconds)
steps=round(3.5/h);  % run for 3.5 seconds

for i=1:steps
    % Graphic of the pendulum
    xpend=[0  L*sin(theta)];
    ypend=[0 -L*cos(theta)];
    subplot(1,2,1)
    plot(xpend,ypend,'o-')
    
    time(i)=i*h;
    thetaRec(i)=theta;
    
    
    axis equal
    axis([-L L -L L])
    subplot(1,2,2)
    plot(time,thetaRec,'r-')
    axis equal
    axis([0 3.5 -5 5])
    pause(h)
    
    alpha=-g/L*sin(theta);     % acceleration term
  
    % Velocity-Verlet Method
    theta=theta + h*omega + 0.5*h*h*alpha;
    alphaNext=-g/L*sin(theta);
    omega=omega + 0.5*h*(alpha+alphaNext);
end 

