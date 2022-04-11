clear;

h=0.1;
x=-4:h:4;
m=length(x);

D=2*eye(m);
D=D-diag(ones(m-1,1),+1);
D=D-diag(ones(m-1,1),-1);
D=D/(2*h^2);

potential=100*ones(m,1);
potential(abs(x)<1)=0;
V=diag(potential);

A=D+V;

[vmat,emat]=eig(A);

K=diag(emat);
disp(K(1:11));

% plot(x,vmat(:,1),'b:',x,vmat(:,2),'g--',x,vmat(:,3),'r-.', ...
%     x,vmat(:,4),'c-',x,vmat(:,5),'m--',x,vmat(:,10),'k-')
% legend('\Psi_1','\Psi_2','\Psi_3','\Psi_4','\Psi_5')

% bar(eig(A))

n=1:5;
fineMesh=1:h:5;

plot(n,K(n),'ro')

% for i=1:5
%     K(i)=K(i)/i^2;
% end
% 
% plot(n,K(n),'ro')
% 
% disp(K(1:5))
