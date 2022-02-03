function ydot=ode_optimal_map(y,s,n,dimension,direction,~,~,~)%,lb,ub,writerObj)
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

global bestCost bestDisp bestEff;

% Because of ODE solver, size of y is [10*dimension 1]
y=reshape(y,[10 dimension]); 

jacobfourier = evaluate_jacobian_fourier(y,s,n,dimension,direction);

%% calcuating ydot so that lagrange equation is zero.
    gradfourier = struct();
    % Reshape jacobdisp and jacobstorke for ODE Solver.
    gradfourier.disp = reshape(jacobfourier.disp,[9*dimension 1]);
    gradfourier.stroke = reshape(jacobfourier.stroke,[9*dimension 1]);

    % Calculate the lagrange multiplier at the current fourier coefficient.
    lambda = pinv(gradfourier.disp)*gradfourier.stroke;

    % The ODE is the gradient of the lagrange equation with slight offset
    % of the lagrange multiplier.
    ydot = gradfourier.stroke-(1-0.02)*lambda*gradfourier.disp;

    % The gradient should move regardless of norm of it.
    ydot = -ydot/norm(ydot);
    if isnan(ydot)
        ydot = zeros(size(ydot));
    end
    
    % The frequency should not change.
    ydot = reshape(ydot,[9 dimension]);
    ydot = [ydot; zeros(1,dimension)];
    ydot = reshape(ydot,[10*dimension 1]);
end