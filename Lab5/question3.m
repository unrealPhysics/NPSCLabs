clear;
% hold on;

h=0.25;
x=-4:h:4;
m=length(x);

D=2*eye(m);
D=D-diag(ones(m-1,1),+1);
D=D-diag(ones(m-1,1),-1);
D=D/(2*h^2);

potential=10*ones(m,1);
potential(abs(x)<1)=0;
V=diag(potential);

A=D+V;

[vmat,emat]=eig(A);

psi=ones(m,1);

for i=1:30
    psi=A*psi;
    psi=psi/max(abs(psi));
end
plot(x,psi,'r-',x,vmat(:,end)./max(abs(vmat(:,end))),'b-')

% B=inv(A);
% psi=ones(m,1);
% err=1;
% j=0;
% OMEGA=[1:0.01:1.2];
% 
% for omega=[1:0.01:1.2]
%     j=j+1;
%     psi=ones(m,1);
%     for i=1:10
%         old=psi;
%         psi=psi*(1-omega) + omega*B*psi;
%         psi=psi/max(abs(psi));
%         errors(i)=sum(abs(old./psi-1));
%         err=errors(end);
%     end
%     ERR(j)=err;
%     clear errors;
% %     plot(x,psi,'ro')
% %     pause()
% end
% 
% semilogy(OMEGA,ERR,'ro-')

% bar(eig(A),1)
% bar(eig(B),0.5)

% hold off;

% 
% psi=ones(m,1);
% I=eye(m);
% omega=1.07;
% C=(1-omega)*I+omega*B;

% for j=1:6
%     old=psi;
%     psi=B*psi;
%     psi=psi/max(abs(psi));
%     error(j)=sum(abs(old./psi-1));
% end
% disp(error)
% 
% psi2=ones(m,1);
% 
% for j=1:6
%     old2=psi2;
%     psi2=C*psi2;
%     psi2=psi2/max(abs(psi2));
%     error2(j)=sum(abs(old2./psi2-1));
% end
% disp(error2)
% plot(x,psi,'ro-',x,psi2,'bo')
% 
% hold on;
% bar(eig(C),1)
% bar(eig(B),0.5)
% hold off;

