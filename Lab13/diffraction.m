% clear

N=2;
a=1;
d=1;
const=1;
u=-5:0.01:5;

E=const*a*pi*sinc(a.*u);

for n=1:(N-1)
    for k=1:length(u)
        shifter=exp(-2*pi*i*(n*d)*u(k));
        E(k)=E(k)+const*a*pi*sinc(a*u(k))*shifter;
    end
end

% EPlot=E/E(floor(length(u)/2)+1);
EPlot=E;
% EPlot=real(EPlot);%+imag(EPlot);
figure(1)
plot(u,real(EPlot),'r-',u,imag(EPlot),'b-')