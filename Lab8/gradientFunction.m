function jacob = gradientFunction(paramVec,const)
%OPT.GRADFNHNDL Non-functional analytic gradient method (incompleted)
%   Non-functional



gravGrad(201:300)=const.effectiveMassMat;

differenceVecX=[const.poles(1,1)-paramVec(1) diff(paramVec(1:10)).' const.poles(2,1)-paramVec(10) diff(paramVec(11:20)).' diff(paramVec(21:30)).' diff(paramVec(31:40)).' diff(paramVec(41:50)).' diff(paramVec(51:60)).' diff(paramVec(61:70)).' diff(paramVec(71:80)).' diff(paramVec(81:90)).' const.poles(3,1)-paramVec(91) diff(paramVec(91:100)).' const.poles(4,1)-paramVec(100)];
differenceVecY=[const.poles(1,2)-paramVec(101) diff(paramVec(101:110)).' const.poles(2,2)-paramVec(110) diff(paramVec(111:120)).' diff(paramVec(121:130)).' diff(paramVec(131:140)).' diff(paramVec(141:150)).' diff(paramVec(151:160)).' diff(paramVec(161:170)).' diff(paramVec(171:180)).' diff(paramVec(181:190)).' const.poles(3,2)-paramVec(191) diff(paramVec(191:200)).' const.poles(4,2)-paramVec(200)];
differenceVecZ=[const.poles(1,3)-paramVec(201) diff(paramVec(201:210)).' const.poles(2,3)-paramVec(210) diff(paramVec(211:220)).' diff(paramVec(221:230)).' diff(paramVec(231:240)).' diff(paramVec(241:250)).' diff(paramVec(251:260)).' diff(paramVec(261:270)).' diff(paramVec(271:280)).' diff(paramVec(281:290)).' const.poles(3,3)-paramVec(291) diff(paramVec(291:300)).' const.poles(4,3)-paramVec(300)];

% springGrad(1:100)=const.kMat*((1-(const.springLengthsNative([181,1:10,182,11:90,183,91:100,184])/const.springLengths([181,1:10,182,11:90,183,91:100,184])))*(differenceVecX(1:100))+((1-(const.springLengthsNative([181,1:10,182,11:90,183,91:100,184])/const.springLengths([181,1:10,182,11:90,183,91:100,184])))*(differenceVecX(1:100))) %x
% springGrad(101:200)= %y
% springGrad(201:300)= %z
% jacob=gravGrad+springGrad;

end

