clear

h=0.05;
x=0:h:1;
y=0:h:1;
m=length(x);
omega=2/(1+sinpi(1/m));

epsilon_0=1;
rho=zeros(m);
rho(round(m*0.4),round(m*0.4))=+1/h^2;
rho(round(m*0.6),round(m*0.4))=+1/h^2;
rho(round(m*0.4),round(m*0.6))=+1/h^2;
rho(round(m*0.6),round(m*0.6))=+1/h^2;
% rho(round(m*0.4),round(m*0.6))=1/h^2;
% rho(round(m*0.6),round(m*0.6))=-1/h^2;
% rho(round(m*0.4),round(m*0.5))=1/h^2;
% rho(round(m*0.6),round(m*0.5))=-1/h^2;
% rho(round(m*0.5),round(m*0.4))=1/h^2;
% rho(round(m*0.5),round(m*0.6))=-1/h^2;

phi=ones(m);
% phi(m,:)=1;
%     for i=2:(m-1)
%         for j=2:(m-1)
%             phi(i,j)=0.25;
%         end
%     end


for n=1:100
    figure(1)
    mesh(x,y,phi)
    title(n)
    drawnow

%     newphi=phi;
    for i=2:(m-1)
        for j=2:(m-1)
            phi(i,j)=(1-omega)*phi(i,j)+ ...
                     (omega/4)* ...
                     (phi(i-1,j)+phi(i,j-1) ...
                     +phi(i+1,j)+phi(i,j+1) ...
                     +h^2*rho(i,j)/epsilon_0);
        end
    end
%     phi=newphi;
end

figure(2)
a=min(min(phi));
b=max(max(phi));
c=linspace(a,b,50);
contour(x,y,phi,c)
axis equal

figure(3)
[ex,ey]=gradient(phi,x,y);
quiver(x,y,ex,ey)
axis equal

interPhi=interp(phi(:),5);
interPhi=reshape(interPhi,m,[]);
figure(4)
a=min(min(interPhi));
b=max(max(interPhi));
c=linspace(a,b,50);
contour(interp(x,5),interp(y,5),interPhi,c)
axis equal