% clear
format compact

% Extent of the System
h=0.1;
x=0:h:1;
y=0:h:1;
m=length(x);
% omega=2/(1+sin(pi/m));
omega=1;
epsilon_0=1;
critVal=1.1;

% Initial Conditions 
phi=zeros(m);

rho=zeros(m);
% for i=1:m
%     for j=1:m
%         critPoint=16*((x(i)-0.5)^2+(y(j)-0.5)^2);
%         if critPoint<critVal && critPoint>1/critVal
%             rho(i,j)=1/h^2;
%         end
%     end
% end
% rho(round(m/2):end,round(m/2)-1:round(m/2)+1)=0;

rho(round(m/2),round(m/2))=1/h^2;

% phi=rand(m);
% for i=1:length(x)
%     for j=1:length(j)
%         phi(i,j)=4/(pi*sinh(pi)) * sin(pi*x(i)) * sinh(pi*y(j));
%     end
% end
% phi(1,:)=0;
% phi(:,1)=0;
% phi(:,end)=0;

% Boundary Conditions: phi=0 everywhere except for y=1
% phi(m,:)=1;

% Main loop
for n=1:1e4
  if rem(n,50)==0
    mesh(x,y,phi)
    title(n)
    drawnow
  end

  % Compute average of surrounding mesh points
  next=phi;
  for i=2:m-1
    for j=2:m-1
      phi(i,j)=(1-omega)*phi(i,j)+omega*0.25*(phi(i-1,j)+phi(i,j-1)+ ...
                      phi(i+1,j)+phi(i,j+1)+ ...
                      h^2*rho(i,j)/epsilon_0);
    end
  end
  
  phi(2:(m-1),1)=phi(2:(m-1),2);
  phi(2:(m-1),m)=phi(2:(m-1),m-1);
  phi(1,2:(m-1))=phi(2,2:(m-1));
  phi(m,2:(m-1))=phi(m-1,2:(m-1));
  phi(1,1)=phi(1,2);
  phi(m,1)=phi(m,2);
  phi(1,m)=phi(2,m);
  phi(m,m)=phi(m-1,m);

  % Evaluate convergence criteria
  diff=abs(next-phi);
  if max(max(diff))<1e-6
    break
  end
    
  % Update to the new value of phi
%   phi=next;
end
figure(4)
mesh(x,y,phi)
title(n)

figure(2)
contour(x,y,phi)
% iterNum(end+1)=n;
% N(end+1)=m;