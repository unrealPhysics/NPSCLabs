clear

h=0.1; % timestep of 0.1s
pos=[1,0];
vel=[0,1];
steps=500;

%frameskipper
skip=ceil(steps/100);

% Circular Orbit
theta=linspace(0,2*pi,100);
xc=cos(theta);
yc=sin(theta);

% Initialise arrays for "speed"
x(steps)=0;
y(steps)=0;

for i=1:steps
    x(i)=pos(1);
    y(i)=pos(2);
    time=(i-1)*h/(2*pi);

    t(i)=time;
    total(i)=0.5*norm(vel)^2 - 1/norm(pos);

    if (rem(i,skip)==0)
        subplot(2,1,1)
        plot(xc,yc,'b',x,y,'g-',pos(1),pos(2),'ko',0,0,'ro')
        title(time)
        axis equal
        drawnow

        subplot(2,1,2)
        plot(t,total);
    end

    r=norm(pos);
    accel=-1/r^2 * pos/r;
    
    % Integrate

%     % Midpoint method (too much stacking error)
%     pos=pos + h*vel + 0.5*h^2*accel;
%     vel=vel + h*accel;

%     % Verlet method
%     if (i==1)
%         posPrev=pos;
%         pos=pos + h*vel + 0.5*h^2*accel;
%         vel=vel + h*accel;
%     else
%         pos=2*pos - posPrev + h^2*accel;
%         vel=(pos-posPrev)/(2*h);
%         posPrev=[x(i),y(i)];
%     end
    

    % Velocity Verlet method
    pos=pos + h*vel + 0.5*h^2*accel;
    accelNext=-1/r^2 * pos/r;
    vel=vel + 0.5*h*(accel+accelNext);

end
