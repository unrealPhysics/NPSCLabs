function CoefftMat = BuildTidalLSQCoefftMat(DayNum, PeriodDays)
%Builds the coeeficient (A) matrix for the tidal fitting least squares
%problem
%DayNum is a vector of decimal days corresponding to height measurements
%PeriodDays = is a length NConstit vector containing the periods (in days) of all the tidal constituents
%
%CoefftMat is the resulting coefficient matrix.
%The first column is all 1s.
%The next NConstit columns contain the sin terms.
%The final NConstit columns contain the cosine terms
%This ordering makes it very easy to build the matrix.
%
%Alec Duncan, Curtin University, 6/4/2019

NDat = length(DayNum);
Omega = 2*pi./PeriodDays;

%Build coefficient matrix:
Mat1 = ones(NDat, 1);
[OmegaMat, DayMat] = meshgrid(Omega, DayNum);

Mat2 = sin(OmegaMat .* DayMat);
Mat3 = cos(OmegaMat .* DayMat);

CoefftMat = [Mat1 Mat2 Mat3];