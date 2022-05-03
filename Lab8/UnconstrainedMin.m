function [ParamVec, FnVal, ParamTrace, FnTrace, RetStatus] =  UnconstrainedMin(FnHndl, Param0, Method, ExtraInfo, Opt)
%Function to implement a number of different unconstrained optimisation
%algorithms to minimise a function for demo purposes
%
%Usage:
%[ParamVec, MinFnVal, ParamTrace, FnTrace] =  UnconstrainedMin(FnHndl, Param0, Method, ExtraInfo, Opt)
%
%Inputs:
%FnHndl - Handle to function to be minimised (name of function with an @ in
%  front of it)  This function must return a scalar value to be minimised
%  when called with the syntax:
%  FnVal = FnHndl(ParamVec, ExtraInfo)
%  where ParamVec  is a column vector of parameter values
%
%Param0 - Initial guess of parameter vector
%Method - 'Gradient descent', 'Conjugate gradient', 'Newton', 'BFGS'
%ExtraInfo - Any data type, but usually a structure, that is passed to the function to be minimised
%   to provide it with whatever extra information it needs
%Opt - optional structure containing optional parameters that control the optimisation.
%These are:
%  .MaxIt - Will abort if max iterations exceeds this number [1000]
%  .AbsTol - Will stop when changes in FnVal from one iteration to the next
%      are below this [1e-6]
%  .RelTol - Will stop when relative changes in FnVal from one iteration to the next are below this
%      [1e-4]
%  .GradFnHndl - If this is a handle to a function then it will be called
%      to calculate the gradient (Jacobian) otherwise this will be done
%      numerically.
%
%Returns:
% ParamVec - Vector of parameter values at the minimum or empty if the
%   minimum couldn't be found
% FnVal - function value at the minimum or empty if the minimum couldn't
%   be found
%ParamTrace - NParam x NIt matrix of parameter values for each iteration of
%  the algorithm
%FnTrace - 1 x NIt vector of function values at each iteration
%RetStatus - 0 if the funciton converged on a minimum, 1 if it failed to
%   reach specified tolerance in Opt.MaxIt iterations, 2 if the function
%   used to carry out the steps had an error
%
%Alec Duncan
%12/4/2020

NPar = length(Param0);

%Sort out defaults for the control parameters.  Doing it this way means the
%user only has to provide any fields in the Opt structure they want to make
%different from the defaults
if nargin < 5
    Opt = [];
end

%Will abort if max iterations exceeds this
if ~isfield(Opt, 'MaxIt')
    Opt.MaxIt = 1000;
end

%Will stop when changes in FnVal from one iteration to the next are below this
if ~isfield(Opt, 'AbsTol')
    Opt.AbsTol = 1e-6;
end

%Will stop when relative changes in FnVal from one iteration to the next are below this
if ~isfield(Opt, 'RelTol')
    Opt.RelTol = 1e-4;
end

if ~isfield(Opt, 'GradFnHndl')
    Opt.GradFnHndl = [];
end

%Initialisations
Converged = false;
Abort = false;
ItCount = 1;
ParamVec = Param0;
FnVal = FnHndl(Param0, ExtraInfo);


%Make space for the trace information
if nargout > 2
    ParamTrace = zeros(NPar, Opt.MaxIt+1);
    FnTrace = zeros(1, Opt.MaxIt+1);
    ParamTrace(:, 1) = ParamVec;
    FnTrace(1) = FnVal;
end

    switch Method
        case 'Gradient descent'
            StepFnHndl = @lSteepestDescentStep;
            
        case 'Conjugate gradient'
            StepFnHndl = @lConjGradientStep;
            
        case 'Newton'
            StepFnHndl = @lNewtonStep;
            
        case 'BFGS'
            StepFnHndl = @lBFGSStep;
    end
    
    StepFnHndl('Reset');  %Reset any internal data
    
%Main loop - the basic mechanics of this are the same for all methods
while ~Converged && ~Abort
    [NewParamVec, NewFnVal] = StepFnHndl(ParamVec, FnVal, FnHndl,  ExtraInfo, Opt);
    
    if isempty(NewParamVec)
        Abort = true;
        RetStatus = 2;
        warning('Step functionr eturned empty parameter vector');
    end
    
    if ~Abort
        %Test for convergence etc.
        dFn = abs(NewFnVal - FnVal);
        if (dFn <= Opt.AbsTol) || (dFn/abs(FnVal) <= Opt.RelTol)
            Converged = true;
            RetStatus = 0;
        elseif ItCount >= Opt.MaxIt
            Abort = true;
            RetStatus = 1;
            warning('Max iterations exceeded');
        end
        %Do updates for the next iteration
        %Not strictly necessary if Aborting but this way I can include the
        %final values in the trace data
        ParamVec = NewParamVec;
        FnVal = NewFnVal;
        ItCount = ItCount + 1;
        if nargout > 2
            ParamTrace(:, ItCount) = ParamVec;
            FnTrace(ItCount) = FnVal;
        end
        bodies=reshape(ParamVec,[],3);
    
        plot3(bodies(:,1),bodies(:,2),bodies(:,3),'k.-',ExtraInfo.poles(:,1),ExtraInfo.poles(:,2),ExtraInfo.poles(:,3),'rx')
        drawnow
    end
end

if nargout > 2
    %Delete any unused memory from the trace arrays
    if ItCount < Opt.MaxIt
        ParamTrace(:, ItCount+1:end) = [];
        FnTrace(ItCount+1:end) = [];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [NewParamVec, NewFnVal] = lSteepestDescentStep(ParamVec, FnVal, FnHndl, ExtraInfo, Opt)
%Steepest descent (gradient descent) method

if ischar(ParamVec) && strcmpi(ParamVec, 'Reset')
    %No initialisations required
else
    %Vector defining search direction is in direction of negative gradient
    if isempty(Opt.GradFnHndl)
        DirnVec = -lNumericalJacobian(ParamVec, FnHndl, ExtraInfo).';
    else
        DirnVec = -Opt.GradFnHndl(ParamVec, ExtraInfo);
    end
    
    %Now we need to search for the minimum along this direction
    %In order to use fminbnd for this we have to set up an anonymous function
    %handle (Yuck!).
    LineMinFnHndl = @(Lambda) FnHndl(ParamVec+Lambda*DirnVec, ExtraInfo);
    
    %We know the lower limit on Lambda is zero, but need to find an upper limit
    LambdaMax = lFindUpperLim(LineMinFnHndl, FnVal);
    if isempty(LambdaMax)
        NewParamVec = [];
        NewFnVal = [];
    else
        BndOpt = optimset('TolX', LambdaMax/100); %Don't want it to try too hard
        [Lambda, NewFnVal] = fminbnd(LineMinFnHndl, 0, LambdaMax, BndOpt);
        NewParamVec = ParamVec+Lambda*DirnVec;
    end
end

val=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [NewParamVec, NewFnVal] = lConjGradientStep(ParamVec, FnVal, FnHndl, ExtraInfo, Opt)
%Conjugate gradient method

%Use persistent variables to carry info over to the next iteration
persistent LastConjDirn
persistent LastSteepestDirn
persistent ResetCount;


if ischar(ParamVec) && strcmpi(ParamVec, 'Reset')
    LastConjDirn = [];
    LastSteepestDirn = [];
    ResetCount = 0;
else
    if ~isfield(Opt, 'NReset')
        Opt.NReset = 20;  %Reset the conjugate direction calculation every this many iterations
    end
    
    %Vector defining search direction is in direction of negative gradient
    if isempty(Opt.GradFnHndl)
        SteepestDirn = -lNumericalJacobian(ParamVec, FnHndl, ExtraInfo).';
    else
        SteepestDirn = -Opt.GradFnHndl(ParamVec, ExtraInfo);
    end
    
    if isempty(LastConjDirn)
        ConjDirn = SteepestDirn;
    else
        %Polak-Ribiere formula
        BetaPR = SteepestDirn.' * (SteepestDirn - LastSteepestDirn) / ...
            (LastSteepestDirn.' * LastSteepestDirn);
        Beta = max(BetaPR, 0);
        
        %Fletcher-Reeves formula
%         Beta = SteepestDirn.' * SteepestDirn / (LastSteepestDirn.' * LastSteepestDirn);
        
        ConjDirn = SteepestDirn + Beta * LastConjDirn;
    end
    
    %Now we need to search for the minimum along this direction
    %In order to use fminbnd for this we have to set up an anonymous function
    %handle (Yuck!).
    LineMinFnHndl = @(Lambda) FnHndl(ParamVec+Lambda*ConjDirn, ExtraInfo);
    
    %We know the lower limit on Lambda is zero, but need to find an upper limit
    LambdaMax = lFindUpperLim(LineMinFnHndl, FnVal);
    if isempty(LambdaMax)
        NewParamVec = [];
        NewFnVal = [];
    else
        BndOpt = optimset('TolX', LambdaMax/100); %Don't want it to try too hard
        [Lambda, NewFnVal] = fminbnd(LineMinFnHndl, 0, LambdaMax, BndOpt);
        NewParamVec = ParamVec+Lambda*ConjDirn;
    end
    LastSteepestDirn = SteepestDirn;
    ResetCount = ResetCount + 1;
    if ResetCount > Opt.NReset
        ResetCount = 0;
        LastConjDirn = [];
    else
        LastConjDirn = ConjDirn;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [NewParamVec, NewFnVal] = lNewtonStep(ParamVec, FnVal, FnHndl, ExtraInfo, Opt)
%Steepest descent (gradient descent) method
    
if ischar(ParamVec) && strcmpi(ParamVec, 'Reset')
    %No initialisations required
else
    if ~isfield(Opt, 'NewtonDamping')
        Opt.NewtonDamping = 1;  %Newton method damping factor, 1 = exact Newton method
    end
    
    if isempty(Opt.GradFnHndl)
        Jac = lNumericalJacobian(ParamVec, FnHndl, ExtraInfo).';
    else
        Jac = Opt.GradFnHndl(ParamVec, ExtraInfo);
    end
    
    Hess = lNumericalHessian(ParamVec, FnHndl, ExtraInfo);
    
    %Hess\Jac is more efficient than inv(Hss) * Jac;
    NewParamVec = ParamVec - Opt.NewtonDamping * Hess\Jac;
    NewFnVal = FnHndl(NewParamVec, ExtraInfo)
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [NewParamVec, NewFnVal] = lBFGSStep(ParamVec, FnVal, FnHndl, ExtraInfo, Opt)
%Steepest descent (gradient descent) method

persistent InvHessian
persistent Grad;
if ischar(ParamVec) && strcmpi(ParamVec, 'Reset')
    InvHessian = [];
    Grad = [];
else
    NPar = length(ParamVec);
    if isempty(InvHessian)
        InvHessian = eye(NPar);
    end
    
    if isempty(Grad)
        if isempty(Opt.GradFnHndl)
            LastGrad = lNumericalJacobian(ParamVec, FnHndl, ExtraInfo).';
        else
            LastGrad = Opt.GradFnHndl(ParamVec, ExtraInfo);
        end
        %Vector defining search direction is in direction of negative gradient
    else
        LastGrad = Grad;
    end
    DirnVec = -InvHessian * LastGrad;

    %Now we need to search for the minimum along this direction
    %In order to use fminbnd for this we have to set up an anonymous function
    %handle (Yuck!).
    LineMinFnHndl = @(Lambda) FnHndl(ParamVec+Lambda*DirnVec, ExtraInfo);
    
    %We know the lower limit on Lambda is zero, but need to find an upper limit
    LambdaMax = lFindUpperLim(LineMinFnHndl, FnVal);
    if isempty(LambdaMax)
        NewParamVec = [];
        NewFnVal = [];
    else
        BndOpt = optimset('TolX', LambdaMax/100); %Don't want it to try too hard
        [Lambda, NewFnVal] = fminbnd(LineMinFnHndl, 0, LambdaMax, BndOpt);
        Sk = Lambda * DirnVec;
        NewParamVec = ParamVec+Sk;
        
        if isempty(Opt.GradFnHndl)
            Grad = lNumericalJacobian(NewParamVec, FnHndl, ExtraInfo).';
        else
            Grad = Opt.GradFnHndl(NewParamVec, ExtraInfo);
        end
        
        Yk = Grad - LastGrad;
        
        %Update the estimate of the Hessian
        SkTYk = Sk.' * Yk;
        
        InvHessian = InvHessian + ...
            (SkTYk + Yk.'*InvHessian*Yk)*(Sk*Sk.') / (SkTYk^2) - ...
            (InvHessian*Yk*Sk.' + Sk*Yk.'*InvHessian)/SkTYk;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Grad = lNumericalJacobian(Params, FnHndl, ExtraInfo)

%Computes the Jacobian of the specified function of Params numerically
NPar = length(Params);

Minda = 1e-10;
da = abs(Params)*1e-8;
IsSmall = da < Minda;
da(IsSmall) = Minda;
LowerParams = Params - da;
UpperParams = Params + da;

%ThisVec = FnHndl(Params, ExtraInfo);

for Ind = 1:NPar
    ThisUpper = Params;
    ThisLower = Params;
    
    ThisUpper(Ind) = UpperParams(Ind);
    ThisLower(Ind) = LowerParams(Ind);
    
    D = (FnHndl(ThisUpper, ExtraInfo) - ...
        FnHndl(ThisLower, ExtraInfo))/ (2*da(Ind));
    
    if Ind == 1
        Grad = zeros(length(D), NPar);
    end
    Grad(:, Ind) = D;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Hess = lNumericalHessian(Params, FnHndl, ExtraInfo)

%Computes the Hessian of the specified function of Params numerically
%By differentiating the Jacobian
NPar = length(Params);

Minda = 1e-10;
da = abs(Params)*1e-4;
IsSmall = da < Minda;
da(IsSmall) = Minda;
LowerParams = Params - da;

Jac0 = lNumericalJacobian(Params, FnHndl, ExtraInfo);

Hess = zeros(NPar, NPar);
for Ind = 1:NPar
    ThisParams = Params;
    ThisParams(Ind) = LowerParams(Ind);
    JacI = lNumericalJacobian(ThisParams, FnHndl, ExtraInfo);
    Hess(Ind, :) = (Jac0 - JacI)/da(Ind);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LambdaMax = lFindUpperLim(LineMinFnHndl, FnVal)
%Aim is to find a point that we know is on the other side of the minimum
%When we know Lambda has to be positive.  FnVal is the function value at Lambda = 0
%This is the dodgiest bit of the whole process
Lambda = 0.01;
Factor = 2;
LastVal = FnVal;
Done = false;
Abort = false;
MaxIt = 1000;
ItCount = 1;
while ~Done && ~Abort
    NewVal = LineMinFnHndl(Lambda);
    if NewVal > LastVal
        LambdaMax = Lambda;
        Done = true;
    else
        ItCount = ItCount + 1;
        
        if ItCount > MaxIt
            Abort = true;
        else
            Lambda = Lambda * Factor;
            LastVal = NewVal;
        end
    end
end

if Abort
    LambdaMax = [];
    warning('Failed to find upper bound');
end

