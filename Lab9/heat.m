clear

kappa=1;
tau=1e-4;
h=0.02;
steps=10;

x=0:h:1;
y=0:h:1;
m=length(x);

D=-2*eye(m);
D=D+diag(ones(m-1,1),-1);
D=D+diag(ones(m-1,1),+1);
% e=ones(m,1);
% D=spdiags([e,-2*e,e],-1:1,m,m);
speedParam=kappa*tau/h^2; % unstable at and above 0.5
disp(speedParam)
D=D*speedParam;
clear e;

D(1,1,:)=0;
D(1,m,:)=0;
D(m,1,:)=0;
D(m,m,:)=0;
temp=zeros(m,m);
temp(round(m/2),round(m/2))=1/h;

time=0;

for n=1:steps
    mesh(x,y,temp)
    title(time)
    drawnow

    time=time+tau;
    temp=temp+pagemtimes(temp;

end