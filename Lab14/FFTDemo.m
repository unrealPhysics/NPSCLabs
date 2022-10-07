function FFTDemo()
%Simple program to demonstrate the use of Matlab's FFT function
%Alec Duncan, Physics and Astronomy, Curtin University
%August 2019

%Step 1 - Simulate a time domain signal to work with.  This is going to
%consist of two sinusoids in noise.
%Define the signal parameters:
TMax = 0.3;  %Signal duration (s)
F1 = 100; %Frequency of first sinusoid
A1 = 1.0; %Amplitude of first sinusoid

F2 = 120; %Frequency of second sinusoid
A2 = 0.5; %Amplitude of second sinusoid
A2=0

dT = 1e-3;  %Sample interval (s)
Fs = 1/dT;  %Sample rate (Samples/s)

% NoiseRms = 1.0;  %RMS noise
NoiseRms=0;

TVec = 0:dT:TMax;  %Build time vector
NPt = length(TVec);  %Number of points in signal

CleanSig = A1*sin(2*pi*F1*TVec) + A2*sin(2*pi*F2*TVec); %Noise free signal
Sig = CleanSig + NoiseRms * randn(1, NPt);  %Add the noise

%Now zero pad to a length equal to the next power of 2.
%Calculating the next power of two is easy if you realise that 
%log_a(X) = log_b(X)/log_b(a) (where a and b are the logarithm bases)
%I use this to calculate log_2(NPt) and then use ceil() to round up to the
%next integer
NextPower = ceil(log10(NPt)/log10(2));
NFFT = 2^NextPower; %Desired FFT length
NFFT

%Matlab's fft() function can zero-pad the data out to any desired length for you but I'm going
%to do it here so I can plot the resulting time series
NPad = NFFT - NPt; %Number of zeros required
PaddedSig = [Sig zeros(1,NPad)];

%I'm also going to calculate a padded time vector purely for plotting
TPad = TVec(end) + (1:NPad) * dT;
PaddedTVec = [TVec TPad];

%Plot the noisy signal, clean signal, and padded signal
figure(1);
clf;
hold on;
plot(TVec, Sig, 'b');
plot(PaddedTVec, PaddedSig, 'r--');
plot(TVec, CleanSig, 'k', 'LineWidth', 1);
xlabel('Time (s)')
ylabel('Signal');
legend('Noisy',  'Zero-padded', 'Clean');
grid on;
box on;

%Now do the FFT - this is the easy bit
Spec = fft(PaddedSig);

%Now calculate the corresponding frequency vector. Note that this has to be
%based on the padded signal length, NFFT, which is all that fft knows about, 
%not the original signal length (NPt).  
%There are different ways to do this, but this one works for me
dF = Fs/NFFT;  %Frequency interval
FPos = (0:NFFT/2)*dF;  %Positive frequencies
%Negative frequency vector excludes first point (0 frequency) and last
%point (NFFT/2*dF). fliplr (flip left-right) reverses the order of the elements in the
%vector.
FNeg = -fliplr(FPos(2:end-1)); 
FVec = [FPos FNeg];

%Plot the result.  Bear in mind that Spec will be complex
%Doing it this way makes sense if you are going to do processing in the frequency domain and then use ifft
%to transform back to the time domain but when plottting you get annoying "fly-back" lines
%between the max and min frequencies because of the order of the data
figure(2);
lPlotSpectrum(FVec, Spec, 'Original order');

%You can avoid the fly-back lines by re-ordering the vectors so the negative 
%frequencies come first.  Matlab's fftshift function is supposed to do
%this, but it swaps the first and second halves of a vector, which is not 
%correct as you actually need to
%put the last NFFT/2-1 points before the first NFFT/2+1 points.
FPlt = [FNeg FPos];
SpecPlt = [Spec(NFFT/2+2:end) Spec(1:NFFT/2+1)];

figure(3);
lPlotSpectrum(FPlt, SpecPlt, 'Re-ordered');

%If you are only interested in the positive frequencies you can just extract those:
IndPos = 1:NFFT/2+1;  
SpecPos = Spec(IndPos);

figure(4);
lPlotSpectrum(FPos, SpecPos, 'Positive frequencies');

%Spectra often have a wide dynamic range, so 
%it is often useful to plot the magnitude in dB and use a logarithmic
%frequency scale
figure(5);
lPlotdBSpectrum(FPos, SpecPos, 'dB spectrum')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function lPlotSpectrum(FVec, Spec, TitleStr)

clf;
subplot(2,1,1);
hold on;
plot(FVec, real(Spec));
plot(FVec, imag(Spec));
plot(FVec, abs(Spec));
legend('Real', 'Imaginary', 'Maginitude');
ylabel('Spectrum');
grid on;
box on;
title(TitleStr);


subplot(2,1,2);
plot(FVec, rad2deg(angle(Spec)));
ylabel('Phase (deg)');
xlabel('Frequency (Hz)');
grid on;
box on;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function lPlotdBSpectrum(FVec, Spec, TitleStr)
clf;
hold on;
plot(FVec, 20*log10(abs(Spec)));
ylabel('Spectrum (dB)');
set(gca, 'XScale', 'log');
grid on;
box on;
title(TitleStr);

xlabel('Frequency (Hz)');

