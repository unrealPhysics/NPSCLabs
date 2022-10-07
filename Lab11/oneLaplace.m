% clear

% Define grid
h=0.1;
x=0:h:1;

tolerance=0.001;
omega=1.2;
% omega=1.55;
% omega=2;

% Determine matrix dimensions
m=length(x);

phi=zeros(m,1);
phi(1)=1;
phi(m)=2;
% phi=phi;

A=zeros(m);
A(1,1)=1;
A(m,m)=1;
for n=2:m-1
    A(n,:)=0.5*omega*A(n-1,:);
    A(n,n)=A(n,n)+(1-omega);
    A(n,n+1)=A(n,n+1)+0.5*omega;
end

iterCount=0;
while sum(abs(phi.'-x-1))/m>tolerance
    iterCount=iterCount+1;
    phi=A*phi;
end

plot(x,phi,'ro-')
disp(iterCount)
% iterLog(end+1)=iterCount;
% omegaLog(end+1)=omega;

sort(abs(eig(A^20)))
sort(abs(eig(A)))