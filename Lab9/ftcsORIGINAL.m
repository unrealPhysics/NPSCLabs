clear

kappa=1;   % thermal conductivity
tau=1e-4;  % time-step
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
while (time<=0.002)
    plot(x,temp,'ro-')   % plot the temperature
    title(time)          % display current time
    pause(0.5)

    time=time+tau;   % update the time
    temp=A*temp;     % update the temperature
end
