clear

h=0.2;
x=-4:h:4;
m=length(x);

D=2*eye(m);
D=D-diag(ones(m-1,1),+1);
D=D-diag(ones(m-1,1),-1);
D=D/(2*h^2);

V=0.5*diag(x.^2);

A=D+V;


