clear;

% h=0.2;
% x=-4:h:4;
% m=length(x);

k=1;
% hold off
% 
% newplot
% 
% hold on

i=0;

hvect=[0.01,0.1,0.25,0.5,0.8,1.0];
logvect=log(hvect);

for h=hvect
    x=-4:h:4;
    m=length(x);
    i=i+1;

%     D=30*eye(m);
%     D=D+diag(ones(m-2,1),+2);
%     D=D-diag(16*ones(m-1,1),+1);
%     D=D-diag(16*ones(m-1,1),-1);
%     D=D+diag(ones(m-2,1),-2);
%     D=D/(12*h^2);
D=2*eye(m);
D=D-diag(ones(m-1,1),+1);
D=D-diag(ones(m-1,1),-1);
D=D/(2*h^2);
    
    V=0.5*k*diag(x.^2);
    
    A=D+V;
    
    [vmat,emat]=eig(A);

    error(i)=abs(emat(1)-0.5);
    herror(i)=error(i)*(h^2);
    logerror(i)=log(error(i));

    
% 
%     if k==1
%         plot(x,vmat(:,1),'r-.');
%     elseif k==2
%         plot(x,vmat(:,1),'g--');
%     elseif k==3
%         plot(x,vmat(:,1),'b:');
%     else
%         plot(x,vmat(:,1),'c-');
%     end
% 
%     sum=0;
%     for i=1:m-1
%         sum=sum+0.5*h*(vmat(i)+vmat(i+1));
%     end
%     disp(sum)
%     title(k)
%     pause(2)
% 
%     if k>1
%         fraction(k-1)=groundeigenvalues(k)/groundeigenvalues(1)
%     end
% 
%     if k>1
%         invfrac(k-1)=1/fraction(k-1)
%     end

end

% plot(hvect,error,'ro-')
% plot(logvect,logerror,'ro-')
plot(hvect.^2,error,'ro-')
for j=1:i-1
    gradient(j)=(logerror(j)-logerror(j+1))/(logvect(j)-logvect(j+1))
end


% for k=1:3
%     deviance(k)=groundeigenvalues(k+1)-groundeigenvalues(k)*sqrt(1+1/k)
%     pattern(k)=((deviance(k)-groundeigenvalues(k+1))/groundeigenvalues(k))^2-1
% end
% 
% legend('k=1','k=2','k=3','k=4')