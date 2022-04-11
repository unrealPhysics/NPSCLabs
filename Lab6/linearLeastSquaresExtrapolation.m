clear;

i=0;

for constituentCount = [5,37]
    i=i+1;
    filepath='./H051011A_2018.csv';
    refYear=2018;
    
    [dayNum,tideHeight]=ReadQueenslandTideData(filepath,refYear);
    
    estVar=(0.001)^2;
    
    weightMat=speye(length(dayNum))/estVar;
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

    goodnessOfFit=std(residualVec);
    
%     subplot(2,1,1)
%     plot(dayNum,tideHeight,'r.-',dayNum,calcVec,'b.')
%     title(chiSquare)
% %     axis([0,10,0,4])
    
%     subplot(2,1,2)
%     plot(dayNum,residualVec,'r.')
%     title(unitVariance)
% %     axis([0,2,-1,1])
%     pause()
    
    % Second Year

    filepath2='./H051011A_2019.csv';
    
    [dayNum2,tideHeight2]=ReadQueenslandTideData(filepath2,refYear);
    
    coeffMat2=BuildTidalLSQCoefftMat(dayNum2,periodDays);
    calcVec2=coeffMat2 * thetaVec;
    residualVec2=calcVec2-tideHeight2;
    
    figure(i)

    subplot(2,1,1)
    plot(dayNum,tideHeight,'r.',dayNum,calcVec,'b-', ...
        dayNum2,tideHeight2,'r.',dayNum2,calcVec2,'g-')
    axis([0,731,0,4])

    title(sprintf('The standard deviation of 2018 residuals at %d constituent waves is %f\n',constituentCount,goodnessOfFit))

    subplot(2,1,2)
    plot(dayNum,residualVec,'b.',dayNum2,residualVec2,'g.')
    axis([0,731,-0.7,0.7])

    goodnessOfFit2=std(residualVec2);
    title(sprintf('The standard deviation of 2019 residuals at %d constituent waves is %f\n',constituentCount,goodnessOfFit2))

end




% n=length(dayNum);
% m=length(periodDays);
% 
% yMat=tideHeight';
% AMat=zeros(2m+1,n); % 2m+1 x n matrix
% Mat1=ones(1,n);
% Mat2=ones(n,m);
% Mat3=ones(n,m);