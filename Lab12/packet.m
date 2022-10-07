% clear

% Define mesh in space and time
h=0.005;
tau=0.25e-4;
x=h:h:1;
m=length(x);
steps=round(0.005/tau);

% Second derivative matrix
D=-2*eye(m);
D=D+diag(ones(m-1,1),+1);
D=D+diag(ones(m-1,1),-1);

% Periodic Boundary Conditions
D(1,m)=1;
D(m,1)=1;

% Define the matrix A - note that by default i=sqrt(-1) 
A=0.5*i*tau/h^2 * D;

% Crank-Nicolson method
mat=inv(eye(m)-A/2) * (eye(m)+A/2);

% Implicit FTCS
% mat=inv(eye(m) - A);

% Explicit FTCS
% mat=eye(m)+A;

% Initial conditions of the wavefunction
% k=100;        % average wavenumber (average wavelength=2*pi/k)
sigma=0.05    % width-parameter for the Gaussian

psi=exp(i*k*x');                           % oscillatory component
psi=psi.* exp(-0.5*((x-0.25)/sigma)'.^2);  % Gaussian envelope
psi=psi/sqrt(sigma*sqrt(pi));              % normalisation

% Remember maximum magnitude for sensible plotting
maxpsi=max(abs(psi));

prob=conj(psi).*psi;
startPos=x(max(prob)==prob);

for n=1:steps
    % Plot real & imaginary components of the wave-packet
%     subplot(2,1,1)
%     plot(x,real(psi),'r',x,imag(psi),'b')
%     axis([0 1 -maxpsi maxpsi])
%     title(n*tau)
% 
%     % Plot probability distribution of the wave-packet
%     subplot(2,1,2)
    prob=conj(psi).*psi;
%     plot(x,prob,'o-')
%     axis([0 1 0 maxpsi^2])
%     drawnow

%     P_tot(n)=h*sum(prob);
%     time(n)=n*tau;
    
    % Advance to next step
    psi=mat*psi;
end

pos=x(max(prob)==prob);
time=n*tau;

