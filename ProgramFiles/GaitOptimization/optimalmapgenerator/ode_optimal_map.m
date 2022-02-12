function ydot=ode_optimal_map(y,s,n,dimension,direction,nfparam)%,lb,ub,writerObj)
%%%%%%%%%%%%%
% This function calculates efficiency (or displacement, if
% that is the objective function) and its gradient with respect to the coefficients obtained
% by the fourier series parametrization
% 
% Inputs: 
% 
% y: Matrix containing the Fourier series coefficients
% s: System file which contains the connection vector field, CCF's and
%   metric data
% dimension: Indicates the number of shape variables of the system
% n: The number of points desired in a direct transcription parametrization
%   of the gaits
% direction: direction in which to optimize motion: 1-x, 2-y, 3-theta
% costfunction: costfunction type to optimize for
% lb: Lower bound of shape variables for each point which is obtained from the grid inside 
%   which an optimal gait is desired
% ub: Upper bound of shape variables for each point which is obtained from the grid inside 
%   which an optimal gait is desired
% 
% Outputs:
% 
% ydot :
%%%%%%%%%%%%%
 
% Because of ODE solver, size of y is [10*dimension 1]
y=reshape(y,[nfparam dimension]); 

[jacobfourier,disp,cost] = evaluate_jacobian_fourier(y,s,n,dimension,direction);

%% Record the fourier coefficient at step-optimal gaits.

global currentDisp bestDisp stepOptimalGaits;

curDispPer= disp/bestDisp;
prevDispPer = currentDisp/bestDisp;
if abs(curDispPer - 3/4) < abs(prevDispPer - 3/4)...
    && abs(curDispPer - 3/4) < 1e-2
    stepOptimalGaits{2,1} = reshape(y,[nfparam dimension]);
    stepOptimalGaits{2,2} = [disp,cost];
    stepOptimalGaits{2,3} = [jacobfourier.disp jacobfourier.stroke];
elseif abs(curDispPer - 2/4) < abs(prevDispPer - 2/4)...
    && abs(curDispPer - 2/4) < 1e-2
    stepOptimalGaits{3,1} = reshape(y,[nfparam dimension]);
    stepOptimalGaits{3,2} = [disp,cost];
    stepOptimalGaits{3,3} = [jacobfourier.disp jacobfourier.stroke];
elseif abs(curDispPer - 1/4) < abs(prevDispPer - 1/4)...
    && abs(curDispPer - 1/4) < 1e-2
    stepOptimalGaits{4,1} = reshape(y,[nfparam dimension]);
    stepOptimalGaits{4,2} = [disp,cost];
    stepOptimalGaits{4,3} = [jacobfourier.disp jacobfourier.stroke];
end

currentDisp = disp;

%% calcuating ydot so that lagrange equation is zero.
gradfourier = struct();

% Reshape jacobdisp and jacobstorke for ODE Solver.
for i = 1:dimension
    gradfourier.disp = reshape(jacobfourier.disp,[(nfparam-1)*dimension 1]);
    gradfourier.stroke = reshape(jacobfourier.stroke,[(nfparam-1)*dimension 1]);
end
% Calculate the lagrange multiplier at the current fourier coefficient.
lambda = pinv(gradfourier.disp)*gradfourier.stroke;

% The ODE is the gradient of the lagrange equation with slight offset
% of the lagrange multiplier.
ydot = gradfourier.stroke-(1-0.02)*lambda*gradfourier.disp;

% The gradient should move toward shirinking the parameter.
ydot = -ydot;
if isnan(ydot)
    ydot = zeros(size(ydot));
end

% The frequency should not change.
ydot = reshape(ydot,[nfparam - 1 dimension]);
ydot = [ydot; zeros(1,dimension)];
ydot = reshape(ydot,[nfparam*dimension 1]);
end