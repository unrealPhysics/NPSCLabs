function FilterSim()
%Example for NPSC2001 lecture that simulates the effect of an electronic
%band-pass filter.
%
%Alec Duncan, 
%Physics and Astronomy
%Curtin University
%27/8/2019

ZXconj=0;
XXconj=0;

SigTypes = {'White noise', 'Pulse', 'Sweep'};
FilterTypes = {'RC lowpass',  'Brick wall lowpass', 'Raised cosine lowpass', ...
    'RLC bandpass', 'Brick wall bandpass', 'Raised cosine bandpass' };  

dT = 1e-5;   %Simulation time step (s)
TSig = 0.1;  %Signal duration (s)
NSig = round(TSig/dT); %Number of samples
SpecDynRange = 60;  %Dynamic range for spectrum plots (dB)

%Make sure we have an even number of points
if rem(NSig, 2) ~= 0
    NSig = NSig + 1;
end

NFFT = NSig;
SigmaSignal = 1;   %Signal rms amplitude


%Choose filter type
FilterOpt = menu('Filter type:', FilterTypes);
FilterType = FilterTypes{FilterOpt};
   
%Choose signal type
SigOpt = menu('Input signal type:', SigTypes);
SigType = SigTypes{SigOpt};
    
mexperiment=20;
for experimentNum=1:mexperiment
%Create the required signal
TVec = (0:NSig-1)*dT;
switch SigType
    case 'White noise' 
        XSig = randn(1, NSig) * SigmaSignal;  %Input signal
    case 'Pulse'
        PulseDelay = 5e-3; %(s)
        PulseLength = 10e-3; %(s)
        XSig = zeros(1, NSig);
        InSig = (TVec >= PulseDelay) & (TVec <= PulseDelay+PulseLength);
        XSig(InSig) = SigmaSignal;
    case 'Sweep'
        %Instantaneous angular frequency is Omega = dPhi/dt
        %So the argument of the sin is Phi = Integral(Omega).dt
        %This results in a division by 2 in the argument that you might not expect
        FMin = 100;
        FMax = 40000;
        Arg = 2*pi*(FMin*TVec + (FMax - FMin)*TVec.^2 / (2*TVec(end)));
        XSig = sqrt(2)*SigmaSignal * sin(Arg);
end
    

%Calculate the filter transfer function for positive frequencies
%Do this at the frequencies we need for the FFT
dF = 1/(NSig*dT);
FVec1 = (0:NSig/2)*dF;         %Frequency vector for positive frequencies
Omega = 2*pi*FVec1;

switch FilterType    
    case 'RC lowpass'
        %Circuit is a series resistor and a parallel capacitor
        Rs = 200000;   %Series resistance (Ohms)
        Cp = 2e-9;     %Parallel capacitance (Farads)
        Av = 1./(1 + 1i*Omega*Rs*Cp);
        
    case 'Brick wall lowpass'
        Fc = 400;  %Cutoff frequency (Hz)
        Av = zeros(1, length(FVec1));
        %Easiest to do this using logical indexing
        IsBelowCutoff = (FVec1 <= Fc);
        Av(IsBelowCutoff) = 1;
        
    case 'Raised cosine lowpass'
        Fc = 400;  % Actually the - 6 dB point
        Av = zeros(1, length(FVec1));
        UseCosine = FVec1 <= 2*Fc;
        Av(UseCosine) = 0.5*(1 + cos(pi*FVec1(UseCosine)/(2*Fc)));
        
    case 'RLC bandpass'
        %Circuit is a series resistor shunted by a parallel tuned circuit - same as first PHYS2004 lab
        Rs = 200000;   %Series resistance (Ohms)
        Rp = 150000;   %Parallel resistance (Ohms)
        Lp = 0.3;      %Parallel inductance (Henrys)
        Cp = 10e-9;     %Parallel capacitance (Farads)
        Av = 1i*Omega*Lp/Rs ./ ((1 - Omega.^2*Lp*Cp) + 1i*Omega*Lp*(1/Rp + 1/Rs));


    case 'Brick wall bandpass'
        FMin = 2800;
        FMax = 3000;
        Av = zeros(1, length(FVec1));
        %Easiest to do this using logical indexing
        IsInBand = (FVec1 >= FMin) & (FVec1 <= FMax);
        Av(IsInBand) = 1; 
        
    case 'Raised cosine bandpass'
        FMin = 2800;  %These will be the -6 dB points
        FMax = 3000;
        Av = zeros(1, length(FVec1));
        F0 = (FMin + FMax)/2;
        BWidth = FMax - FMin;
        IsInBand = (FVec1 >= F0-BWidth) & (FVec1 <= F0+BWidth);
        Av(IsInBand) = 0.5 *(1 + cos(pi*(FVec1(IsInBand) - F0)/BWidth));
end

%Make the full frequency vector that includes the negative frequencies
%and the corresponding transfer function.  Remember the negative frequency
%values have to be complex conjugates of the corresponding positive
%frequency values.
FVecAll = [FVec1 -fliplr(FVec1(2:end-1))];
AvAll = [Av conj(fliplr(Av(2:end-1)))];


XSpec = fft(XSig);  %Calculate the spectrum of the signal

YSpec = XSpec .* AvAll;  %Multiply by the transfer function

YSig = ifft(YSpec);  %Finally inverse fft to get the fitered signal

YSigma = 0.001; % measurement noise

YSig = YSig + YSigma.*randn(1,length(YSig));

ZSpec = fft(YSig);

TitleTxt = ['Input: ' SigType ', Filter: ' FilterType];

%Plot the filter transfer - first on linear axes
figure(1);
clf;
hold on;
plot(FVec1, abs(Av))
ylabel('Transfer function (linear scale)');
title(TitleTxt);
xlabel('Frequency (Hz)');
grid on;
box on;


figure(2);
clf;
hold on;
PltVals = 20*log10(abs(Av));
MaxPltVal = ceil(max(PltVals)/10)*10; %Round to nxt highest multiple of 10 dB
MinPltVal = MaxPltVal-SpecDynRange;
IsBad = isinf(PltVals) | isnan(PltVals);
PltVals(IsBad) = MinPltVal;
plot(FVec1, PltVals)
ylabel('Transfer function (dB)');
title(TitleTxt);
xlabel('Frequency (Hz)');
set(gca, 'XScale', 'log')
set(gca, 'YLim', [MinPltVal MaxPltVal]);
grid on;
box on;


%Plot up the time domain signals
figure(3);
clf
subplot(2,1,1);
hold on;
plot(TVec, XSig);
grid on;
box on;
ylabel('Input signal');
xlabel('Time (s)');
title(TitleTxt);


subplot(2,1,2);
hold on;
%YSig should theoretically be real but rounding error will result in it
%having a small imaginary part.  It is a good idea to plot this to make
%sure it is small.  If it isn't it means there is a problem with the way
%the transfer function has been extended to include negative frequencies.
%This has saved me many times, including this time!
plot(TVec, real(YSig));
plot(TVec, imag(YSig));
legend('Output - real', 'Output - imaginary');
grid on;
box on;
ylabel('Output signal');
xlabel('Time (s)');




%And their spectra (positive frequencies only)
IndPos = 1:NFFT/2+1;

figure(4);
clf;
hold on;
PltVals = 20*log10(abs(XSpec(IndPos)));
MaxPltVal = ceil(max(PltVals)/10)*10; %Round to next highest multiple of 10 dB
MinPltVal = MaxPltVal-SpecDynRange;

YPltVals = 20*log10(abs(YSpec(IndPos)));
IsBad = isnan(YPltVals) | isinf(YPltVals);
YPltVals(IsBad) = MinPltVal;

plot(FVec1, PltVals)
plot(FVec1, YPltVals)
ylabel('Spectrum (dB)');
title(TitleTxt);
xlabel('Frequency (Hz)');
legend('Input', 'Output')
set(gca, 'XScale', 'log')
set(gca, 'YLim', [MinPltVal MaxPltVal]);
grid on;
box on;

ZXconj=ZXconj+ZSpec.*conj(XSpec);
XXconj=XXconj+XSpec.*conj(XSpec);

end

aveZXconj=ZXconj/mexperiment;
aveXXconj=XXconj/mexperiment;

figure(5)
clf;
hold on;
plot(FVec1,real(aveZXconj(IndPos)))
plot(FVec1,imag(aveZXconj(IndPos)))
plot(FVec1,real(aveXXconj(IndPos)))
plot(FVec1,imag(aveXXconj(IndPos)))
legend('ZX* - real', 'ZX* - imaginary', 'XX* - real', 'XX* - imaginary');
grid on;
box on;

for n=1:length(aveZXconj)
    H(n)=aveZXconj(n)/aveXXconj(n);
end

figure(6)
clf;
hold on;
plot(FVec1,abs(H(IndPos)),'r-')
plot(FVec1,angle(H(IndPos)),'g-.')
plot(FVec1,abs(Av),'c--')
plot(FVec1,angle(Av(IndPos)),'m:')
ylabel('Transfer function (linear scale)');
title(TitleTxt);
xlabel('Frequency (Hz)');
legend('H - absolute', 'H - phase', 'Av - absolute', 'Av - phase');
axis([min(FVec1),max(FVec1),-1,1])
grid on;
box on;


