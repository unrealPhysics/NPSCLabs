function [value,R,U] = shootingSchrodinger(e)

r=100;
h=0.0001;
steps=r/h;

u=0;
prev=1;

for i=1:steps
  R(i)=r;
  U(i)=u;

  d2u = -2 * (u/r + e*u);
  next=2*u - prev + h*h*d2u;

  prev=u;
  u=next;

  r=r-h;
end

area=h*sum(U.^2);
U=U./sqrt(area);
%PSI=PSI/max(abs(PSI));
plot(R,U)
title(e)
drawnow

value=U(end);

