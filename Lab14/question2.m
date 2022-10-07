clear

N=2;
a=1;
d=10;
magE_0=1;
h=0.02;

% phi=-pi/2:pi/100:pi/2;

y=[zeros(1,d/h) ones(1,a/h)];

for n=2:N
    y=[y zeros(1,d/h) ones(1,a/h)];
end
y=[y zeros(1,d/h)];

x=h:h:N*(a+d)+d;

totalLength=(length(x)-1)*h;
spatialFreq=((1:length(x))-2)/(2*totalLength);

E_0=magE_0*y;%*exp(i*phi);

spec=fft(E_0);

figure(1)
clf
loglog(spatialFreq,abs(spec),'r-')

