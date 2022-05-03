%DEPRECATED DO NOT USE

const.g=9.8;
longBoi=ones(10,1);
const.effectiveMassMat=const.g*ones(100,1);
bodies(:,1)=reshape(longBoi*[1:10],[],1); %#ok<*NBRAK> 
bodies(:,2)=reshape([1:10].'*longBoi.',[],1);
bodies(:,3)=zeros(100,1);
const.poles=[0,0,5;0,11,5;11,0,5;11,11,5];
const.k=0.5;
const.kMat=const.k*0.5*ones(184,1);
const.springLengthsNative=ones(184,1);


potentialGravitational=const.effectiveMassMat.'*bodies(:,3);

for i=1:10
    springLengths(((i-1)*9+1):(9*i))=(diff(bodies((((i-1)*10+1):(10*i)),1)).^2 ...
        +diff(bodies((((i-1)*10+1):(10*i)),2)).^2 ...
        +diff(bodies((((i-1)*10+1):(10*i)),3)).^2).^(0.5);
    springLengths((90+i):10:180)=(diff(bodies(i:10:100,1)).^2 ...
        +diff(bodies(i:10:100,2)).^2 ...
        +diff(bodies(i:10:100,3)).^2).^(0.5);
end

springLengths(181)=sqrt((bodies(1,1)-const.poles(1,1))^2+ ...
                        (bodies(1,2)-const.poles(1,2))^2+ ...
                        (bodies(1,3)-const.poles(1,3))^2);
springLengths(182)=sqrt((bodies(10,1)-const.poles(2,1))^2+ ...
                        (bodies(10,2)-const.poles(2,2))^2+ ...
                        (bodies(10,3)-const.poles(2,3))^2);
springLengths(183)=sqrt((bodies(91,1)-const.poles(3,1))^2+ ...
                        (bodies(91,2)-const.poles(3,2))^2+ ...
                        (bodies(91,3)-const.poles(3,3))^2);
springLengths(184)=sqrt((bodies(100,1)-const.poles(4,1))^2+ ...
                        (bodies(100,2)-const.poles(4,2))^2+ ...
                        (bodies(100,3)-const.poles(4,3))^2);

potentialElastic=((springLengths-const.springLengthsNative.').^2*const.kMat);
totalPotential=potentialElastic+potentialGravitational;


% gravGradX=0;
% gravGradY=0;
% gravGradZ=const.effectiveMassMat;


% [ParamVec, FnVal, ParamTrace, FnTrace, RetStatus] = UnconstrainedMin(FnHndl, Param0, Method, ExtraInfo, Opt)