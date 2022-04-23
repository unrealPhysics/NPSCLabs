% H=newplot;
% reset(H);
clear;
format compact;
% H=newplot;

noiseStandardDev=0.01;

converged=false;
aborted=false;

repLimit=1000;
epsilon=0.00001;

lambda=1; % IF ZERO, SIMPLE METHOD USED. IF POSITIVE, LEVENBURG-MARQUARDT METHOD USED. IF NEGATIVE, PROCESS WILL LIKELY FAIL.
nu=1.5;
chiSquareOld=0;
pointCount=10;

error1=noiseStandardDev*randn(1,pointCount);

upperLimit=10;

truePoint=[34;-12];

knownCoords=upperLimit*rand(2,pointCount)-0.5*[upperLimit;upperLimit];

estimatedPoint=[30;-15];

for noiseStandardDev=[0.1,1,2,5,10,20]
    converged=false;
    aborted=false;
    
    repLimit=1000;
    epsilon=0.00001;
    
    lambda=1; % IF ZERO, SIMPLE METHOD USED. IF POSITIVE, LEVENBURG-MARQUARDT METHOD USED. IF NEGATIVE, PROCESS WILL LIKELY FAIL.
    nu=1.5;
    chiSquareOld=0;
    pointCount=10;
    
    error1=noiseStandardDev*randn(1,pointCount);
    
%     upperLimit=10;
%     
%     truePoint=[34;-12];
%     
%     knownCoords=upperLimit*rand(2,pointCount)-0.5*[upperLimit;upperLimit];
%     
    estimatedPoint=[30;-15];
    % estimatedPoint=truePoint+[truePoint(1)*(rand(1)-0.5);truePoint(2)*(rand(1)-0.5)] % decent estimate mode
    
    WMat=speye(pointCount)/((noiseStandardDev)^2);
    
    d=sqrt((truePoint(1)-knownCoords(1,:)).^2+(truePoint(2)-knownCoords(2,:)).^2);
    dMeas=d+error1;
    
    repCount=0;

    while ~converged && ~aborted
        % Repetition Counting
        repCount=repCount+1;
    
        % Parameter Construction
        r=sqrt((estimatedPoint(1)-knownCoords(1,:)).^2+(estimatedPoint(2)-knownCoords(2,:)).^2);
    
        rMeas=r+error1;
    
        yf=dMeas.'-rMeas.';
    
        for i=1:pointCount
            AMat(i,1)=(estimatedPoint(1)-knownCoords(1,i))/rMeas(i);
            AMat(i,2)=(estimatedPoint(2)-knownCoords(2,i))/rMeas(i);
        end
    
        DMat=diag(AMat.' * WMat * AMat);
        deltaVec=(AMat.' * WMat * AMat + lambda * DMat)\(AMat.' * WMat * yf);
        estimatedPointLambda=estimatedPoint+deltaVec;
    
        if lambda~=0
            r=sqrt((estimatedPointLambda(1)-knownCoords(1,:)).^2+(estimatedPointLambda(2)-knownCoords(2,:)).^2);
            rMeas=r+error1;
            residual=(rMeas-dMeas).';
            chiSquareLambda=(residual.' * WMat * residual);
    
            if chiSquareOld==0
                chiSquareOld=chiSquareLambda;
            end
            
            deltaVecNu=(AMat.' * WMat * AMat + lambda/nu * DMat)\(AMat.' * WMat * yf);
            estimatedPointNu=estimatedPoint+deltaVecNu;
    
            r=sqrt((estimatedPointNu(1)-knownCoords(1,:)).^2+(estimatedPointNu(2)-knownCoords(2,:)).^2);
            rMeas=r+error1;
            residual=(rMeas-dMeas).';
            chiSquareNu=(residual.' * WMat * residual);
    
            if chiSquareLambda > chiSquareOld
                if chiSquareNu < chiSquareOld
                    chiSquareOld=chiSquareNu;
                    lambda=lambda/nu;
                    estimatedPoint=estimatedPointNu;
                    deltaVec=deltaVecNu;
                else
                    lambda=nu*lambda;
                end
            elseif chiSquareNu < chiSquareLambda
                chiSquareOld=chiSquareNu;
                lambda=lambda/nu;
                estimatedPoint=estimatedPointNu;
                deltaVec=deltaVecNu;
            else
                chiSquareOld=chiSquareLambda;
                estimatedPoint=estimatedPointLambda;
            end
        else
            estimatedPoint=estimatedPointLambda;
        end
    
        if norm(deltaVec) < epsilon
            converged=true;
        end
    
        if repCount > repLimit
            aborted=true;
        end
    
        locX(repCount)=estimatedPoint(1);
        locY(repCount)=estimatedPoint(2);
    end
    
    residual=(rMeas-dMeas).';
    chiSquare=(residual.' * WMat * residual);
    unitVar=chiSquare / (length(dMeas)-length(deltaVec));
    covarMat=unitVar*inv(AMat.' * WMat * AMat);
    
    positionString=string(num2str(estimatedPoint(1)));
    positionString(2,1)=string(num2str(estimatedPoint(2)));
    positionString(1,2)=join(['(' num2str(truePoint(1)) ')']);
    positionString(2,2)=join(['(' num2str(truePoint(2)) ')']);
    
    if converged
        disp(['Converged in ' num2str(repCount) ' iterations.'])
        disp('The estimated position (true position) of the unknown point is ')
        disp(positionString)
        disp('The covariance matrix of final iteration is ')
        disp(covarMat)
    elseif aborted
        disp(['Failed to converge in ' num2str(repCount) ' iterations.'])
        disp('The estimated position (true position) of the unknown point is ')
        disp(positionString)
        disp('The covariance matrix of final iteration is ')
        disp(covarMat)
    else
        disp('Severe error in logic! Finished before convergence or abortion!')
    end
end

% for j=1:pointCount
%     changeY=estimatedPoint(2)-knownCoords(2,j);
%     changeX=estimatedPoint(1)-knownCoords(1,j);
%     dist=sqrt(changeX^2+changeY^2);
%     propX=changeX/dist;
%     propY=changeY/dist;
%     distort(2,j,1)=(dMeas(j)*propY+knownCoords(2,j));
%     distort(1,j,1)=(dMeas(j)*propX+knownCoords(1,j));
%     xVal(j)=distort(1,j,1);
%     yVal(j)=distort(2,j,1);
%     knownX(j)=knownCoords(1,j);
%     knownY(j)=knownCoords(2,j);
%     hold on
%     diamondX=[knownX(j)+dMeas(j),knownX(j)+dMeas(j),knownX(j)-dMeas(j),knownX(j)-dMeas(j),knownX(j)+dMeas(j)];
%     diamondY=[knownY(j)+dMeas(j),knownY(j)-dMeas(j),knownY(j)-dMeas(j),knownY(j)+dMeas(j),knownY(j)+dMeas(j)];
%     plot(diamondX,diamondY,'g-')
% end
% % plot(xVal,yVal,'k*',knownCoords(1,j),knownCoords(2,j),'ro')
% aveX=sum(xVal)/length(xVal);
% aveY=sum(yVal)/length(yVal);
% plot(locX,locY,'ko-',truePoint(1),truePoint(2),'cs',estimatedPoint(1),estimatedPoint(2),'k+',xVal,yVal,'k*',aveX,aveY,'b+',knownX,knownY,'ro')
% y=0;
% (x-x(1))^2+(y-y(1))^2=radius^2
% (x-x(2))^2+(y-y(2))^2=radius^2
% (x-x(3))^2+(y-y(3))^2=radius^2
% (x-x(1))^2-(x-x(2))^2=(y-y(2))^2-(y-y(1))^2
% x(1)^2-x(2)^2-2*x*x(1)+2*x*x(2)=y(2)^2-y(1)^2-2*y*y(2)+2*y*y(1)


% if aborted
%     point(1)=((x(1)^2-x(3)^2)*(y(1)-x(2))+(y(1)^2-y(3)^2)*(y(1)-y(2))+(x(2)^2-x(1)^2)*(y(1)-y(3))+(y(2)^2-y(1)^2)*(y(1)-y(3)))/(x(3)-x(1))*(y(1)-y(2))-(x(2)-x(1))*(y(1)-y(3))/(-2);
%     point(2)=(x(1)^2-x(2)^2-2*point(1)*x(1)+2*point(1)*x(2)-y(2)^2+y(1)^2)/(y(2)-y(1))/-2;
% 
%     hold on
%     plot(point(1),point(2),'bd')
% end


% 
% for count=2:repCount-1
%     gradient=(locY(count+1)-locY(count-1))/(locX(count+1)-locX(count-1));
%     gradient=-1/gradient;
%     pointX=[locX(count)-2,locX(count)+2];
%     pointY=[locY(count)-2*gradient,locY(count)+2*gradient];
%     plot(pointX,pointY,'b-')
% end
