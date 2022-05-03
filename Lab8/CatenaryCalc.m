function CatenaryCalc()
%Function to demo minimisation of a function of many parameters.
%
%The problem is as follows:
%Consider N identical masses of mass M, connected together by identical
%springs with unstretched length L0 and spring constant Ks.  The two end
%masses are connected to the top of vertical posts of  heights H1 and H2
%which are separated horizontally by distance XDist.
%
%Write a program that can find the equilibrium positions of the masses when given
%specific values of N, M, L0, Ks, H1 and H2 by minimising the potential
%energy of the system.
%
%Alec Duncan, Curtin University, April 2020

%First save all the parameters that define the problem in a structure so we
%can easily pass it through to the cost function

FigNum = 11;
PlotIntermedSolns = false;

Info.H1 = 4.0;  %Height of first post, (m)
Info.H2 = 2.0;  %Height of second post (m)
Info.XDist = 3.0;  %Horizontal separation between posts (m)
Info.NMass = 100; %Number of masses
Info.Mass = 2.0/Info.NMass;  %Mass of each (kg)

%Calculate the unstretched length as a proportioin of the straight line
%distance between the end points (note there are N+1 springs)
SLDist = sqrt(Info.XDist^2 + (Info.H2 - Info.H1)^2);
Info.L0 = 0.9 * SLDist/(Info.NMass+1);
Info.Ks = 3*Info.NMass;  %Spring constant (N/m)

Info.g = 9.8;  %Acceleration due to gravity (m/s^2)


%Coordinate system is X horizontal, Y positive up

%Set up the vector of initial parameters, this will be all the X values
%followed by all the Y values in one big column vector.  Started with the masses on a straight line
%between the end points.
dX = Info.XDist/(Info.NMass + 1);
X0 = (1:Info.NMass).' * dX;
Y0 = (1:Info.NMass).' * (Info.H2 - Info.H1)/(Info.NMass + 1) + Info.H1;

Par0 = [X0; Y0];

XParMin = X0 - Info.XDist/10;
XParMax = X0 + Info.XDist/10;
YParMin = zeros(size(Y0));
YParMax = Y0 + 0.5;

ParMin = [XParMin; YParMin];
ParMax = [XParMax; YParMax];


NPar = length(Par0);
%Plot this up to make sure we got it right

Opt.NewtonDamping = 1;
Opt.NReset = 500;  %Conjugate gradient is reset every this many iterations
Opt.AbsTol = 1e-7;
Opt.RelTol = 1e-7;
Opt.MaxIt = 20000;
%Handle to function that calculates the gradient of the cost function
%(Jacobian).  Empty to use numerical derivatives
%Opt.GradFnHndl = @lCatenaryGradFn;  
Opt.GradFnHndl = [];

%Method = {'Newton'};
disp(' ');
disp('Method, Seconds taken, Number of iterations, Number of cost fn calls, Number of gradient fn calls, Final energy')
%Method = {'Gradient descent', 'Conjugate gradient', 'Newton', 'BFGS'};
%Method = {'Conjugate gradient', 'BFGS'};
Method = {'Conjugate gradient'};

for IMeth = 1:length(Method)
    figure(FigNum);
    FigNum = FigNum + 1;
    clf;
    subplot(2,1,1);
    hold on;
    %Plot the end points
    plot([0, Info.XDist], [Info.H1, Info.H2], 'ro', 'MarkerFaceColor', 'r');
    lPlotPar(Par0, Info);
    lCatenaryCostFn('Reset');
    if ~isempty(Opt.GradFnHndl)
        Opt.GradFnHndl('Reset');
    end
    tic;
    switch Method{IMeth}
        case 'Reduced gradient'
            Opt.UseHessian = true;
            [ParVec, FnVal, ParamTrace, FnTrace] = LibNDMinimiseRG(Info, [], Par0, @lRGDSCatCostFn, ParMin, ParMax, false, Opt);
        case 'Direction set'
            [ParVec, FnVal, ParamTrace, FnTrace] = LibNDMinimise(Info, [], Par0, @lRGDSCatCostFn, ParMin, ParMax, false, Opt);
            
        case 'Simulated annealing'
            Opt.TFactor = 0.99;
            [ParVec, FnVal, ParamTrace, TempTrace, FnTrace, AllParamTrace, AllFnTrace] = ...
                LibSimulatedAnneal(Info, [], Par0, @lRGDSCatCostFn, XParMin, XParMax, false, Opt);
            %ParamTrace = AllParamTrace;
        otherwise
            
            [ParVec, FnVal, ParamTrace, FnTrace, RetStatus] =  UnconstrainedMin(@lCatenaryCostFn, Par0, Method{IMeth}, Info, Opt);
    end
    ElapsedTime = toc;
    if PlotIntermedSolns
        for ITrace = 1:length(FnTrace)
            lPlotPar(ParamTrace(:, ITrace), Info);
        end
    end
    
    if ~isempty(ParVec)
        plot(ParVec(1:NPar/2), ParVec(NPar/2+1:end), 'ko', 'MarkerFaceColor', 'k')
    end
    
    NIt = length(FnTrace);

    title([Method{IMeth} ', ' num2str(Info.NMass) ' masses']);
    xlabel('X (m)');
    ylabel('Y (m)');
    axis([0 Info.XDist 0 1.1*max(Info.H1, Info.H2)]);
    grid on;
    box on;
    
    subplot(2,1,2);
    plot(FnTrace);
    xlabel('Iteration');
    ylabel('Potential energy (J)');
    %axis([0 inf 20 30]);
    grid on;
    box on;
    drawnow;
    title([num2str(NIt) ' iterations, final energy = ' num2str(FnVal) ' J']);
    CostFnCount = lCatenaryCostFn('GetCount');
    if ~isempty(Opt.GradFnHndl)
        GradFnCount = Opt.GradFnHndl('GetCount');
    else
        GradFnCount = 0;
    end

    disp([Method{IMeth} ', ' num2str(ElapsedTime) ', ' num2str(NIt) ', ' num2str(CostFnCount) ', ' num2str(GradFnCount) ', ' num2str(FnVal)]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function lPlotPar(ParVec, Info)
NPar = length(ParVec);
XPlt = [0; ParVec(1:NPar/2); Info.XDist];
YPlt = [Info.H1; ParVec(NPar/2+1:end); Info.H2];
plot(XPlt, YPlt, '.-');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PE = lCatenaryCostFn(ParVec, Info)
%Stuff so I can count the times this is called
persistent CallCount;
if ischar(ParVec)
    switch ParVec
        case 'Reset'
            CallCount = 0;
            PE = [];
        case 'GetCount'
            PE = CallCount;
    end
else
    %Actual cost function calcs
    CallCount = CallCount + 1;
    
    NPar = length(ParVec);
    XVec = ParVec(1:NPar/2);
    YVec = ParVec(NPar/2+1:end);
    
    %Total gravitational potential energy
    PEGravity = sum(Info.Mass * Info.g * YVec);
    
    dX = [XVec(1); diff(XVec); Info.XDist-XVec(end)];
    dY = [YVec(1)-Info.H1; diff(YVec); Info.H2-YVec(end)];
    
    Dist = sqrt(dX.^2 + dY.^2);
    PESprings = sum(0.5*Info.Ks*(Dist - Info.L0).^2);
    
    PE = PEGravity + PESprings;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function FnVal = lRGDSCatCostFn(ExtraInfo, Y, Params)
%This provides a slightly different interface to the cost funciton for
%compatibility with some pre-existing optimisation routines
FnVal = lCatenaryCostFn(Params, ExtraInfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Grad = lCatenaryGradFn(ParVec, Info)
%Analytic calculation of the gradient of the cost function
%To avoid approximations in the numerical calc and speed things up

%Stuff so I can count the times this is called
persistent CallCount;
if ischar(ParVec)
    switch ParVec
        case 'Reset'
            CallCount = 0;
            Grad = [];
        case 'GetCount'
            Grad = CallCount;
    end
else
    %Actual gradient function calcs
    CallCount = CallCount + 1;
    
    NPar = length(ParVec);
    XVec = ParVec(1:NPar/2);
    YVec = ParVec(NPar/2+1:end);
    
    XAll = [0; XVec; Info.XDist];
    YAll = [Info.H1; YVec; Info.H2];
    
    dX = diff(XAll);
    dY = diff(YAll);
    
    DistAll = sqrt(dX.^2 + dY.^2);
    Fact = 1 - (Info.L0 ./ DistAll);
    
    dUdX = - Info.Ks * (Fact(2:end) .* dX(2:end) - Fact(1:end-1) .* dX(1:end-1));
    
    dUdY = - Info.Ks * (Fact(2:end) .* dY(2:end) - Fact(1:end-1) .* dY(1:end-1)) + ...
        Info.Mass * Info.g;
    
    Grad = [dUdX; dUdY];
end
