clear;

g=9.8;             % acceleration due to gravity
L=0.2482;          % length in metres (gives period of 1 second)

theta=0*pi/180;   % initial angle (in radians)
omega=0;           % initial angular-velocity (in radians/second)
h=0.001;            % timestep (in seconds)
steps=179;  % run for 180 degrees

for i=1:steps
    thetaInitials(i)=i*pi/180;
    boundary=ones(i,1).*1.05;
    theta=thetaInitials(i);
    thetaPrev=theta;
    k=0;
    while(theta<=thetaPrev)
        thetaPrev=theta;
        k=k+1;
        time=k*h;
        period(i)=2*time;

        alpha=-g/L*sin(theta);     % acceleration term
      
        % Velocity-Verlet Method
        theta=theta + h*omega + 0.5*h*h*alpha;
        alphaNext=-g/L*sin(theta);
        omega=omega + 0.5*h*(alpha+alphaNext);
        if(theta>0)
            thetaPrev=theta+h;
        end
    end
    plot(thetaInitials,period,'r-',thetaInitials,boundary,'b-')
    drawnow
end
