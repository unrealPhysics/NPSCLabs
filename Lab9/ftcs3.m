clear

kappa=1/2;   % thermal conductivity
tau=1e-4/2;  % time-step
h=0.02;    % spatial-step

% Extent of the system
x=0:h:1;
m=length(x);

% Matrix for the second-derivative term
D=-2*eye(m);
D=D+diag(ones(m-1,1),+1);
D=D+diag(ones(m-1,1),-1);
D=D*kappa*tau/(h*h);

% Boundary conditions (Dirichlet; i.e. constant value)
D(1,:)=0;
D(m,:)=0;

% The update matrix
A=eye(m) + D;

% Initial conditions
temp=zeros(m,1);
temp(round(m/2))=1/h;

time=0;
% midTemp=temp(round(m/2));
% midTime=time;
% midATemp=(1/sqrt(4*pi*kappa*time));
while (time<=0.02)
    plot(x,temp,'ro-')   % plot the temperature
    title(time)          % display current time
    pause(0)

    time=time+tau;   % update the time
    temp=A*temp;     % update the temperature
%     midTemp(end+1)=temp(round(m/2));
%     midATemp(end+1)=(1/sqrt(4*pi*kappa*time));
%     midTime(end+1)=time;
end

% semilogy(midTime,midTemp,'r-',midTime,midATemp,'b-')
% 
% plot(x,temp,'ro-',x,sin(pi.*x)*tau,'bx--')
% max(temp)