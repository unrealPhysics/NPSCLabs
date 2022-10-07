% clear;

% Parameters for the system
c=1;
tau=0.01;
h=0.01;
steps=ceil(1/(c*tau));     % pulse loops through system once

% Extent of the system
x=h:h:1;       % no mesh point at x=0 due to periodic boundary conditions
m=length(x);

% Diagonal matrix to implement the spatial first-derivative
D=diag(ones(m-1,1),+1);
D=D-diag(ones(m-1,1),-1);

% Additional elements for periodic boundary conditions
D(m,1)=1;
D(1,m)=-1;

% Construct the spatial-averaging matrix
av=0.5*abs(D);

% Construct the matrix for the spatial second-derivative
D2=abs(D)-2*eye(m);

% FTCS update matrix
% mat=eye(m) + c*tau/(2*h)*D;

% Implicit FTCS update matrix
% mat=inv(eye(m)-c*tau/(2*h)*D);

% Crank-Nicolson update matrix
D=c*tau/(2*h)*D;
mat=inv(eye(m)-D/2)*(eye(m)+D/2);

% Lax update matrix
% mat=av + c*tau/(2*h)*D;

% Lax-Wendroff update matrix
% mat=eye(m) + c*tau/(2*h)*D + 0.5*(c*tau/h)^2*D2;

% Initial conditions (Gaussian pulse at x=0.5)
wave=exp(-100*(x-0.5).^2)';
wave0=wave; 

% Main loop
for n=1:steps
    wave=mat*wave;    % Update current values

%     plot(x,wave0,'ro-',x,wave,'bo-')
%     title(n*tau)        
% %     area(n)=h*sum(wave);
% %     time(n)=n*tau;
%     drawnow
end

% max(abs(eig(mat)))
% width=h*sum(wave>0.5*max(wave));

% area=h*sum(wave)
error=sum(h*abs(wave-wave0));