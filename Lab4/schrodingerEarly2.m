clear;

%Boundaries and such
x=-4;
initialx=x;
h=0.001;
steps=-2*x/h;
skip=steps/100;

%Values
hbar=1;
m=1;
k=1;
e=0.81973926;

%Initialisers
psi=1;
psiPrev=1;

% X(steps)=0;
% PSI(steps)=0;

for i=1:steps

    %Plotting
    X(i)=x;
    PSI(i)=psi;

%     if rem(i,skip)==0
%         plot(X,PSI,'ro-');
%         %axis([initialx -initialx 0 1])
%         drawnow
%     end

    % Calculate potential
    if abs(x)<1
        pot=0;
    else
        pot=10;
    end
    d2psi=(2*m/hbar^2)*(pot-e)*psi;

    % Verlet Method
    psiNext=2*psi-psiPrev+h^2*d2psi;

    psiPrev=psi;
    psi=psiNext;
    x=x+h;
end

% Plot with exact solution
% plot(X,PSI,'ro',X,exp(-0.5.*X.^2),'b-')
% legend('Numerical','Exact')
plot(X,PSI,'ro-');
pause(1)
area=h*sum(PSI.^2);
PSI=PSI./sqrt(area);
plot(X,PSI,'ro-');