clear;

hbar=1;
m=1;
k=1;
A=1;
x=1;
E=2.5;

omega=sqrt(k/m);

V=0.5*k*x.^2;

psi=A*(2*x^2-1)*exp(-(x^2*m*omega)/(2*hbar));

doubleDerivPsi=(0.5*k*x^2*psi-E*psi)/(hbar^2/(2*m))