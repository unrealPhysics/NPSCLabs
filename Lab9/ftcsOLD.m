clear

kappa=1;   % thermal conductivity
% kappaList=0:0.5:4;
% for i=1:length(kappaList)
% kappa=kappaList(i);
tau=1e-4;  % time-step
% for tau=[0.5e-4,1e-4,1.5e-4,1.75e-4,2e-4,2.25e-4,2.5e-4]
h=0.02;    % spatial-step

% Extent of the system
x=0:h:1;
% aX=0:h/10:1;
m=length(x);

% Matrix for the second-derivative term
D=-2*eye(m);
D=D+diag(ones(m-1,1),+1);
D=D+diag(ones(m-1,1),-1);
speedParam=kappa*tau/h/h;
D=D*speedParam;

% Boundary conditions (Dirichlet; i.e. constant value)
D(1,:)=zeros(1,m);
D(m,:)=5*ones(1,m);

% The update matrix
A=eye(m) + D;

% Initial conditions
temp=zeros(m,1);
temp(round(m/2))=1/h;
% aTemp=zeros(length(aX),1);
% aTemp(round(length(aX/2)))=1/h;

time=0;
while (time<=(10000*tau))
    plot(x,temp,'ro-')   % plot the temperature
    title(time)          % display current time
%     if rem(floor(time/tau),100)==0
        drawnow
%     end

    time=time+tau;   % update the time
    temp=A*temp;     % update the temperature
%     aTemp=(1/sqrt(4*pi*kappa*time))*exp(-(aX-0.5).^2/(4*kappa*time));
end
% pause()
% eigenvals=eig(A);
% maxEig(i)=max(eigenvals);
% minEig(i)=min(eigenvals);
% domEig(i)=max(abs(eigenvals));
% stabVal(i)=speedParam;
% end
% plot(stabVal,minEig,'bo:',stabVal,maxEig,'rx-',stabVal,domEig,'k+--')
% legend('Minimum Eigenvalue','Maximum Eigenvalue','Magnitude of Dominant Eigenvalue','Location','northwest')