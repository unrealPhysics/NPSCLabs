clear;

[e1,x1,psi1]=shooting(0.5);
[e2,x2,psi2]=shooting(1.5);
[e3,x3,psi3]=shooting(2.5);

plot(x1,psi1,'r-',x2,psi2,'g--',x3,psi3,'b-.')
legend('\Psi_1','\Psi_2','\Psi_3')

% psi4=psi1.*(2.*x1.^2.-1);
% C=max(abs(psi3))/max(abs(psi4))
% psi5=C.*psi4;
% plot(x1,psi4,'r-.',x3,psi3,'b-',x1,psi5,'g--')

p1p2=0.01*sum(psi1.*psi2)
p2p3=0.01*sum(psi2.*psi3)
p3p1=0.01*sum(psi3.*psi1)
p1p1=0.01*sum(psi1.*psi1)
p2p2=0.01*sum(psi2.*psi2)
p3p3=0.01*sum(psi3.*psi3)