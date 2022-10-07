clear
inject=[2,4,8];
for n=1:3
TMax=inject(n);
F0=100;
Fs=1e3;
T=1;
TVec=0:1/Fs:T;
% TMax=4;
signal=[cos(2*pi*F0*TVec) zeros(1,ceil((TMax-T)*Fs))];
TMaxVec=0:1/Fs:TMax;
sigLength=length(signal);
titleString=sprintf(['Frequency of pulse: %G Hz, Sample rate: %G Hz\n' ...
    'Duration of pulse: %G s, Total duration: %G s'],F0,Fs,T,TMax);

% figure(1)
% plot(TMaxVec,signal)

spectra=fft(signal);

FPos=(0:sigLength/2)*Fs/sigLength;
FNeg=-fliplr(FPos(2:end-1));
FTot=[FPos FNeg];

%generate analytic solution
analytic=T/2*(sinc(T*(F0-FPos))+sinc(T*(F0+FPos)));

% extract positives
posIndices=1:sigLength/2+1;
spectraPos=spectra(posIndices);

% display
figure(n)
clf
semilogx(FPos,20*log10(abs(spectraPos)/max(abs(spectraPos))),'r-', ...
         FPos,20*log10(abs(analytic)/max(abs(analytic))),'b-')
title(titleString)
xlabel('Frequency (Hz)')
ylabel('Normalised Spectrum (dB)')
legend('Spectra','Theoretical',Location='best')
axis auto

end
