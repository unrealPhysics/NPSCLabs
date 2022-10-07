function DrawnApertureDiffraction()
%This function demostrates using an image to define an aperture and then
%calculating its diffraction pattern using the 2D FFT
%Alec Duncan, 
%Physics and Astronomy, Curtin University
%29/8/2019

ImFile = 'triangleSlitPart1.png';

FigNum = 111;  %Starting figure


PlotType = 'dB';  %Type of colour scale to use.  'linear' or 'dB'

%For convenience the physical dimensions will be specified in wavelengths of the
%incident light

%Image width and height in wavelengths - this will be used to convert pixel
%coordinates to physical units in wavelengths of the inciodent light
ScreenWidth = 150;
ScreenHeight = 100;

ImDat = imread(ImFile); %Read it

% ImDat = 255*ones(size(ImDat),'uint8')-ImDat;

%Convert to grayscale by converting to double and then summing over the colour channels (planes of the matrix) 
%I've flipped it vertically because images are defined with (0,0) in the top
%left corner
GScaleImage = flipud(sum(double(ImDat), 3));  
GScaleImage = GScaleImage+flipud(GScaleImage);

[NY, NX] = size(GScaleImage);
%Trim edges if required to make size even
if rem(NY, 2) ~= 0
    GScaleImage = GScaleImage(1:end-1, :);
    NY = NY - 1;
end

if rem(NX, 2) ~= 0
    GScaleImage = GScaleImage(:, 1:end-1);
    NX = NX - 1;
end


%Make vectors of aperture coordinates (again in wavelengths)
dX = ScreenWidth/(NX-1);
XVec = (0:NX-1)*dX - ScreenWidth/2;

dY = ScreenHeight/(NY-1);
YVec = (0:NY-1)*dY - ScreenHeight/2;


%Make the corresponding spatial frequency vectors.  We have to do this
%twice so it is worth writing a little helper function
UVec = lMakeSpatialFreqVec(dX, NX);
VVec = lMakeSpatialFreqVec(dY, NY);


%2D FFT of full aperture function
%Build a single 2D aperture function
%Note that rows come first in Matlab indexing, so order of indices is Y then X
%Again, logical indexing makes this easy.

% GScaleImage=GScaleImage-mean(GScaleImage(:));
DPat = fft2(GScaleImage);

%Plot the 2D aperture function
FigNum = FigNum + 1;
figure(FigNum);
clf;
colormap('gray');
pcolor(XVec, YVec, GScaleImage);
shading interp;
xlabel('x/\lambda');
ylabel('y/\lambda');
colorbar;
axis equal;

FigNum = FigNum + 1;
figure(FigNum);
lPlotPattern(UVec, VVec, DPat, '2D FFT', PlotType);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function UVec = lMakeSpatialFreqVec(dX, NX)
%Make up a spatial frequency vector in the usual FFT output order
%Assumes NX is even
dU = 1/(NX*dX);
UPos = (0:NX/2)*dU;
UVec = [UPos -fliplr(UPos(2:end-1))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function lPlotPattern(UVec, VVec, DPat, Description, PlotType)
%Plots the 2D diffraction pattern.  This is a bit messy so definitely worth
%its own function

%The main issue is that everything is still in FFT output order, and if we
%just plot it directly that will mess things up. In theory we could use
%Matlab's fftshift function to reorder things, but it doesn't quite do what
%it is supposed to, so we'll do it ourselves

NU = length(UVec);
%This index vector will re-order things correctly in the X (column)
%dimensio
IndU = [NU/2+2:NU 1:NU/2+1];

NV = length(VVec);
IndV = [NV/2+2:NV 1:NV/2+1];

switch PlotType
    case 'linear'
        PltVals = abs(DPat)/max(abs(DPat(:))); %Normalise so max value is 1
        CRange = [0 1];
        CLabel = '';

    case 'dB'
        PltVals = 20*log10(abs(DPat)/max(abs(DPat(:)))); %Amplitude only and normalise so max value is 1
        CRange = [-60 0];
        CLabel = '(dB)';
end

clf;
colormap('default');
pcolor(UVec(IndU), VVec(IndV), PltVals(IndV, IndU));
caxis(CRange);
shading interp;
xlabel('u\lambda = sin(\theta)');
ylabel('v\lambda = sin(\phi)');
title(Description);
Hnd = colorbar;
Hnd.Label.String = CLabel;
axis([-1 1 -1 1]);  %Just plot the region corresponding to real angles





