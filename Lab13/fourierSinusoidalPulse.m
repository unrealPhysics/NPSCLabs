clear

f_0=1;
T=1;
freq=0:0.01:16;
figure(1)

for f_0=0.1:0.01:8

    for n=1:length(freq)
        fourierGH(n)=sin(pi*T*(f_0-freq(n)))/(2*pi*(f_0-freq(n)))+...
                     sin(pi*T*(f_0+freq(n)))/(2*pi*(f_0+freq(n)));
    end

    plot(freq,fourierGH,'r-')
    axis([min(freq),max(freq),-0.5,1])
    title(f_0)
    pause(0)
end
