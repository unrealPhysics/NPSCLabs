clear;

filepath='./H051011A_2018.csv';
refYear=2018;

[dayNum,tideHeight]=ReadQueenslandTideData(filepath,refYear);

estVar=(0.001)^2;

weightMat=speye(length(dayNum))/estVar;

for constituentCount = [1:1:37]
    periodDays=GetTidalConstituentPeriods(constituentCount);
    coeffMat=BuildTidalLSQCoefftMat(dayNum,periodDays);
    
    thetaVec=(coeffMat' * weightMat * coeffMat)\(coeffMat' * weightMat * tideHeight);
    
    calcVec=coeffMat * thetaVec;
    
    residualVec=calcVec-tideHeight;

    chiSquare=residualVec' * weightMat * residualVec;

    popMean=mean(tideHeight);
    unitVariance=chiSquare/(length(dayNum)-constituentCount);

    covarMat=unitVariance*inv(coeffMat' * weightMat * coeffMat);
    standDev=sqrt(diag(covarMat));

    goodnessOfFit(constituentCount)=std(residualVec);
    
    subplot(2,1,1)
    plot(dayNum,tideHeight,'r.-',dayNum,calcVec,'b.')
    title(sprintf('chi^2 =  %f',chiSquare))
%     axis([0,10,0,4])
    
    subplot(2,1,2)
    plot(dayNum,residualVec,'r.')
    title(sprintf('Unit Variance = %f',unitVariance))
%     axis([0,2,-1,1])
    pause()
end

% n=length(dayNum);
% m=length(periodDays);
% 
% yMat=tideHeight';
% AMat=zeros(2m+1,n); % 2m+1 x n matrix
% Mat1=ones(1,n);
% Mat2=ones(n,m);
% Mat3=ones(n,m);