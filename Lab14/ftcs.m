clear

kazarel=0;

kappa=1;   % thermal conductivity
tau=1e-4;  % time-step
h=0.02;    % spatial-step

% Extent of the system
x=0:h:1;
m=length(x);

% Matrix for the second-derivative term
eigVals=-2*eye(m);
eigVals=eigVals+diag(ones(m-1,1),+1);
eigVals=eigVals+diag(ones(m-1,1),-1);
eigVals=eigVals*kappa*tau/(h*h);

% Boundary conditions (Dirichlet; i.e. constant value)
eigVals(1,:)=0;
eigVals(m,:)=0;

% The update matrix
A=eye(m) + eigVals;

% Initial conditions
% temp=ones(m,1);
% temp(round(m/2))=1/h;
% % temp(10)=1/h;
k=0.1;
temp=(k*x).';

NPeriod=5;

partialTemp=[temp(1:m-1), 2*temp(m)-temp(m:2)];
totalTemp=[partialTemp, partialTemp, partialTemp, partialTemp];

[eigVecs,eigVals]=eig(A,'vector');
[sortedEigVals,index]=sort(abs(eigVals),'descend');
sortedEigVecs=eigVecs(:,index);
figure(1)
clf

subplot(2,3,1)
plot(x,sortedEigVecs(:,1),'ro-', ...
     x,sortedEigVecs(:,2),'go-', ...
     x,sortedEigVecs(:,3),'bo-', ...
     x,sortedEigVecs(:,4),'co-', ...
     x,sortedEigVecs(:,5),'mo-', ...
     x,sortedEigVecs(:,6),'ko-')

totalLength=(length(x)-1)*h;
spatialFreq=((1:length(x))-2)/(2*totalLength);

time=0;
while (time<=20*0.002)
    subplot(2,3,2)
    plot(x,temp,'ro-')   % plot the temperature
    title(time)          % display current time
    pause(0)

    a=sortedEigVecs\temp;

    subplot(2,3,3)
    plot(x,a,'ro-',x,abs(a),'bo-')

    subplot(2,3,4)
    fourTemp=fft(temp);
    plot(spatialFreq,abs(fourTemp),'ro-')

    subplot(2,3,5)
    partialTemp=[temp(1:1:m-1);2*temp(m)-temp(m:-1:2)];
    totalTemp=partialTemp;
    if NPeriod~=1
        for n=2:NPeriod
            totalTemp=[totalTemp; partialTemp];
        end
    end
    spec=fft(totalTemp,(2*NPeriod*(m-1)));
    plot(spatialFreq,abs(spec(1:length(spatialFreq))),'ro-')

    subplot(2,3,6)
%     plot(1:length(totalTemp),totalTemp,'bo-')

    runSum=0;

%     spec=spec/max(abs(spec));
%     fourTemp=fourTemp/max(abs(fourTemp));
    for n=1:length(spatialFreq)
%         runSum=runSum+real(fourTemp(n))*sin(2*pi*spatialFreq(n)*x)+imag(fourTemp(n))*cos(2*pi*spatialFreq(n)*x);
        runSum=runSum+real(spec(n))*sin(2*pi*spatialFreq(n)*x)+imag(spec(n))*cos(2*pi*spatialFreq(n)*x);
    end

    plot(x,runSum,'ro-')

    time=time+tau;   % update the time
    temp=A*temp;     % update the temperature
end


% figure(2)
% plot(x,runSum,'ro-')