function [value,X,PSI] = shootingmodified(e)

x=-4;
h=0.01;
steps=-2*x/h;

psi=1;
prev=0;

for i=1:steps
  X(i)=x;
  PSI(i)=psi;

  pot=0.1*x^4;
  d2psi=2*psi*(pot-e);
  next=2*psi - prev + h*h*d2psi;

  prev=psi;
  psi=next;

  x=x+h;
end

area=h*sum(PSI.^2);
PSI=PSI./sqrt(area);
%PSI=PSI/max(abs(PSI));
plot(X,PSI)
title(e)
drawnow

value=PSI(end);

