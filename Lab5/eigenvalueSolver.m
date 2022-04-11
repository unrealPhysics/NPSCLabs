clear;

h=0.01;
approxInf=4;
x=-approxInf:h:approxInf;
m=length(x);

D=(1/(2*h^2)).*(-diag(ones(m-1,1),1)-diag(ones(m-1,1),-1)+diag(2.*ones(m,1)));
V=0.5.*diag(x.^2);

A=D+V;
E=0.5;

% [vmat,emat]=eig(A);
% 
% plot(x,vmat(:,1),'r-',x,vmat(:,2),'g-',x,vmat(:,3),'b-');

B=inv(A);

psi=zeros(m,1);
psi(1)=1;

omega=1.05;

for i=1:10
    plot(x,psi,'ro-');
    title(i)
    old=psi;
    psi=(1-omega)*psi + omega*B*psi;
    pause(0.5)
    maxVal=max(abs(psi));
    psi=psi./maxVal;

    errors(i)=sum(abs(old./psi-1));
end

plot(x,psi,'ro-');
title(mean(B*psi./psi))

semilogy(errors,'ro-');
title(errors(end))