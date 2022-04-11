clear

h=0.001;         % time-step
x=[1 0 0 0.8];  % initial conditions: [x y vx vy]

steps=20000;             % number of steps
skip=ceil(steps/200);   % display 200 frames

for i=1:steps
  % Fourth-order Runge-Kutta
  F1=deriv(x);
  F2=deriv(x+h/2*F1);
  F3=deriv(x+h/2*F2);
  F4=deriv(x+h*F3);

  x=x+h/6*(F1+2*F2+2*F3+F4);
    
  % Save and plot trajectory
  X(i)=x(1);
  Y(i)=x(2);

  KE(i)=0.5*((x(3))^2+(x(4))^2);
  time(i)=i*h;
  if(i>1)
      runningAverage(i)=(runningAverage(i-1)*(i-1)+KE(i))/i;
  else
      runningAverage(i)=KE(i);
  end
   
  if rem(i,skip)==0
%     subplot(1,2,1)
%     plot(X,Y,'g-',X(end),Y(end),'ko',0,0,'ro')
%     axis equal
%     subplot(1,2,2)
    plot(time,KE,'r-',time,runningAverage,'b--');
    drawnow
  end
end

% -------------------

function value = deriv(vec)

k=0.02;
pos=vec(1:2);
vel=vec(3:4);

r=norm(pos);
accel=-1/r^2 * pos/r -0.01*norm(vel)*vel;

value(1:2)=vel;
value(3:4)=accel;

end
