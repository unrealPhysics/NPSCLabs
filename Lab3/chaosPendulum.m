clear
format compact

global g L
g=9.8;
L=0.2482;

theta=45*pi/180;
x=[theta 0];
kx=[theta 0];
h=0.02;

% % RK Method
% F1=derivOld(x);
% F2=derivOld(x+h*F1);
% x=x+h/2*(F1+F2);

for i=1:100
    % RK4 Method
    F1=derivOld(x);
    F2=derivOld(x+h/2*F1);
    F3=derivOld(x+h/2*F2);
    F4=derivOld(x+h*F3);
    x=x+h/6*(F1+2*F2+2*F3+F4);

    kF1=derivOld(kx);
    kF2=derivOld(kx+h/2*kF1);
    kF3=derivOld(kx+h/2*kF2);
    kF4=derivOld(kx+h*kF3);
    kx=kx+h/6*(kF1+2*kF2+2*kF3+kF4);

    % Graphic of the Pendulum
    xpend=[0 L*sin(x(1))];
    ypend=[0 -L*cos(x(1))];
    kxpend=xpend+[0 L*sin(kx(1))];
    kypend=ypend+[0 -L*cos(kx(1))];
    xvals(i)=xpend(2);
    yvals(i)=ypend(2);
    kxvals(i)=kxpend(2);
    kyvals(i)=kypend(2);
    plot(xpend,ypend,'o-',kxpend,kypend,'ro-')
    axis equal
    axis([-2*L 2*L -2*L 0])
    drawnow
end

disp('x=')
disp(x)
disp('1/e=')
disp(exp(-1))
