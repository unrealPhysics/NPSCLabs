clear;

h=0.05;       % timestep of 0.05 time-units
pos=[1 0];    % initial position
vel=[0 1.1];  % initial velocity

steps=200;              % number of steps
skip=ceil(steps/100);   % plot only 100 frames

for i=1:steps
  x(i)=pos(1);
  y(i)=pos(2);
  time=(i-1)*h/(2*pi);

  % Draw the trajectory and ellipse every skip'th step
  if rem(i,skip)==0
    a=0.5*(max(x)-min(x));   % compute width of ellipse
    b=0.5*(max(y)-min(y));   % compute height of ellipse

    x0=0.5*(max(x)+min(x));  % compute centre (x0,y0) of ellipse
    y0=0.5*(max(y)+min(y));

    % Two kinds of ellipse (tall/thin, or squat/wide)
    if a>b
      e=sqrt(1-b^2/a^2);
      xfocus(1)=x0+sqrt(a^2-b^2);  % foci of squat ellipse
      xfocus(2)=x0-sqrt(a^2-b^2);
      yfocus=[y0 y0];  
    else
      e=sqrt(1-a^2/b^2);
      xfocus=[x0 x0];              % foci of tall ellipse
      yfocus(1)=y0+sqrt(b^2-a^2);
      yfocus(2)=y0-sqrt(b^2-a^2);
    end

    % Show the ellipse every 10 degrees
    angle=linspace(0,2*pi,37);   
    xe=x0+a*cos(angle);
    ye=y0+b*sin(angle);

    plot(pos(1),pos(2),'go',x,y,'g-',xe,ye,'ko',0,0,'ro',xfocus,yfocus,'bo')
    title(sprintf('time=%.2f years  eccentricity=%.3f',time,e))
    axis equal
    pause(0.1)
  end

  r=norm(pos);
  accel=-1/r^(-2) * pos/r;
 
  % Verlet Method
  if i==1
    next=pos + h*vel + 0.5*h*h*accel;
  else
    next=2*pos - prev + h*h*accel;
  end
 
  prev=pos;
  pos=next;
end

