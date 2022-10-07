clear

kappa=1*2;   % thermal conductivity
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
D(1,1)=-1;
D(1,2)=1;
D(m,:)=0;
D(m,m-1)=1;
D(m,m)=-1;

% The update matrix
A=eye(m) + D;

% Initial conditions
temp=zeros(m,1);
temp(round(m/2))=1/h;
% temp(26)=1/h;

time=0;
% midTemp=temp(round(m/2));
% midTime=time;
% midATemp=(1/sqrt(4*pi*kappa*time));
while (time<=(1000*tau))
%     plot(x,temp,'ro-')   % plot the temperature
%     title(time)          % display current time
%     pause(0)

    time=time+tau;   % update the time
%     if rem(time,tau)~=0
%         break
%     end
    temp=A*temp;     % update the temperature
%     temp(1)=temp(2);
%     midTemp(end+1)=temp(round(m/2));
%     midATemp(end+1)=(1/sqrt(4*pi*kappa*time));
%     midTime(end+1)=time;
% disp(temp(25)-temp(27))
% temp=round(temp,14);
end

%     plot(x,temp,'ro-',x,(max(temp)+min(temp))/2-(max(temp)-min(temp))/2*cos(2*pi.*x),'bx--',x,min(temp)+(max(temp)-min(temp))*exp(-16*(x-0.5).^2),'go-')   % plot the temperature
plot(x,temp,'bx-')    
title(time)          % display current time
drawnow

% disp(kappa*tau/h^2)
% semilogy(midTime,midTemp,'r-',midTime,midATemp,'b-')