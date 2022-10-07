% clear
if exist('y','var') ~=  1
    [y, Fs]=audioread("twoMinutes.m4a");
%     y=(y(:,1)+y(:,2))/2;
    y=y(1:end-Fs);
end

% startTime=1;
% endTime=4;
% y=y(1+Fs*startTime:Fs*endTime);

if exist('t','var') ~=  1
    t=(0:length(y)-1)/Fs;
    t=t.';
end

% t=Fs*(25:34);
% sig=y(t(1):t(end));
% tEff=t(1)/Fs:1/Fs:t(end)/Fs;
NFFT=Fs/4;

figure(1)
clf
plot(t,y)

figure(2)
clf
spectrogram(y,hanning(NFFT),[],NFFT,Fs,'yaxis')
% spectrogram(y,'yaxis')
ax=gca;
% ax.YScale='log';
% colormap bone
% ax.ColorScale='log';

figure(3)
clf
pwelch(y,hanning(NFFT),[],NFFT,Fs)
ax=gca;
% ax.XScale='log';

figure(4)
clf
pwelch(y([9*Fs:floor(31*Fs) 63*Fs:floor(76*Fs) 100*Fs:floor(115*Fs)]),hanning(NFFT),[],NFFT,Fs)
ax=gca;
% ax.XScale='log';

% player=audioplayer(y,Fs,16);