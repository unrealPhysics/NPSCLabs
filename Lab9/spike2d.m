clear;

kappa=1;   % thermal conductivity
tau=6.25e-4;  % time-step
h=0.05;    % spatial-step

% Extent of the system
x=0:h:1;
y=0:h:1;
m=length(x);

% Initial conditions
temp=zeros(m,m);
temp(round(m/2),round(m/2))=1/h^2;

D=-4*eye(m^2);
D(1,2)=1;

for i=1:m^2
    if rem(i,m)~=0
        if i~=m^2
            D(i,i+1)=1;
        end
    end
    if rem(i,m)~=1
        if i~=1
            D(i,i-1)=1;
        end
    end
    if i-m>0
        D(i,i-m)=1;
    end
    if i+m<m^2
        D(i,i+m)=1;
    end
    imagesc(D)
    pause(0)
end

D=D*kappa*tau/h/h;
imagesc(D)
pause(0.5)

% Boundary Condition
for i=1:m
    D(i,:)=0;
    D(m^2+1-i,:)=0;
    D(m*(i-1)+1,:)=0;
    D(m*i,:)=0;
    imagesc(D)
    pause(0)
end

imagesc(D)
pause(1)
A=eye(m^2)+D;
imagesc(A)
pause(1)

tempVec=temp(:);

time=0;
for n=1:100
    mesh(x,y,temp)      % mesh plot
    title(n*tau);       % display current time
    drawnow

    tempVec=A*tempVec;
    temp=reshape(tempVec,m,[]);

%     next=temp;
%     for i=2:m-1
%       for j=2:m-1
%         factor= temp(i-1,j)+temp(i+1,j)+temp(i,j-1)+temp(i,j+1) -4*temp(i,j);
%         next(i,j)=temp(i,j) + kappa*tau/h^2 * factor;
%       end
%     end
%     temp=next;

end

sort(eig(A))
disp(h^2/4/kappa-tau)