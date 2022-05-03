function Lab8FUNCTIONVERSION()
%LAB8FUNCTIONVERSION Summary of this function goes here
%   Detailed explanation goes here
const.g=9.8;
longBoi=ones(10,1);
const.effectiveMassMat=const.g*ones(100,1);
const.poles=[-1,-1,10;-1,12,10;12,-1,10;12,12,10];
% const.poles=[-1,-1,10;12,-1,10;-1,12,10;12,12,10]; %flipped mode
const.k=1000;
const.kMat=const.k*0.5*ones(184,1);
% const.kMat(1:11:180)=5*const.kMat(1:11:180); %modified spring constants
% const.kMat(1:3:180)=const.kMat(1:3:180)/5; %modified spring constants
const.springLengthsNative=ones(184,1);
% const.effectiveMassMat(55)=-5*const.effectiveMassMat(55); %modified mass

method={'Gradient descent', 'Conjugate gradient', 'BFGS'};
for methodVal=1:length(method)
    bodies(:,1)=reshape(longBoi*[1:10],[],1); %#ok<*NBRAK> 
    bodies(:,2)=reshape([1:10].'*longBoi.',[],1);
    bodies(:,3)=10*ones(100,1);
    % Opt.GradFnHndl=@gradientFunction; % do not reactivate, function is incomplete
    [ParamVec, FnVal, ParamTrace, FnTrace, RetStatus] = UnconstrainedMin(@generatePotential, bodies(:), string(method(methodVal)), const);
    bodies=reshape(ParamVec,[],3);
    
    plot3(bodies(:,1),bodies(:,2),bodies(:,3),'k.-',const.poles(:,1),const.poles(:,2),const.poles(:,3),'rx')
    fprintf('Method be %s, number of iterations be %s, final value be %s\nPress enter to continue\n',string(method(methodVal)),num2str(length(FnTrace)),num2str(FnVal))
    pause()
end

end

function totalPotential = generatePotential(paramVec, const)
bodies = reshape(paramVec,[],3);
potentialGravitational=const.effectiveMassMat.'*bodies(:,3);
% 
% for i=1:10
%     springLengths(((i-1)*9+1):(9*i))=(diff(bodies((((i-1)*10+1):(10*i)),1)).^2 ...
%         +diff(bodies((((i-1)*10+1):(10*i)),2)).^2 ...
%         +diff(bodies((((i-1)*10+1):(10*i)),3)).^2).^(0.5);
%     springLengths((90+i):10:180)=(diff(bodies(i:10:100,1)).^2 ...
%         +diff(bodies(i:10:100,2)).^2 ...
%         +diff(bodies(i:10:100,3)).^2).^(0.5);
% end
    % same as commented above, but non-iterative
springLengths=([1,1,1]*([diff(bodies(1:10,:)).' diff(bodies(11:20,:)).' diff(bodies(21:30,:)).' diff(bodies(31:40,:)).' diff(bodies(41:50,:)).' diff(bodies(51:60,:)).' diff(bodies(61:70,:)).' diff(bodies(71:80,:)).' diff(bodies(81:90,:)).' diff(bodies(91:100,:)).' diff(bodies(1:10:100,:)).' diff(bodies(2:10:100,:)).' diff(bodies(3:10:100,:)).' diff(bodies(4:10:100,:)).' diff(bodies(5:10:100,:)).' diff(bodies(6:10:100,:)).' diff(bodies(7:10:100,:)).' diff(bodies(8:10:100,:)).' diff(bodies(9:10:100,:)).' diff(bodies(10:10:100,:)).'].^2)).^0.5;

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
end