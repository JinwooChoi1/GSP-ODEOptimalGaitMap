function xdot=ode_optimal_map_ellipse_gait(x,s,n,dimension,direction)
%%%%%%%%%%%%%
% This function calculates efficiency (or displacement, if
% that is the objective function) and its gradient with respect to the coefficients obtained
% by the fourier series parametrization
% 
% Inputs: 
% 
% x: r1 r2 for elliptical axis
% s: System file which contains the connection vector field, CCF's and
%   metric data
% dimension: Indicates the number of shape variables of the system
% n: The number of points desired in a direct transcription parametrization
%   of the gaits
% direction: direction in which to optimize motion: 1-x, 2-y, 3-theta
% 
% Outputs:
% 
% xdot :
%%%%%%%%%%%%%

%% Calculating cost and displacement per gait

% Assign a time period for executing the gait
T = 2*pi;

% Define phi_def = [alpha1, alpha2] as a function of time t such that the
% array returns the shape variables given by the fourier coefficients at a
% time t
len_x=length(x);

p = make_ellipse_gait(x);
dx = 1e-6;
if len_x == 2
    dp1 = make_ellipse_gait([x(1)+dx;x(2)]);
    dp2 = make_ellipse_gait([x(1);x(2)+dx]);
elseif len_x == 3
    dp1 = make_ellipse_gait([x(1)+dx;x(2);x(3)]);
    dp2 = make_ellipse_gait([x(1);x(2)+dx;x(3)]);
    dp3 = make_ellipse_gait([x(1);x(2);x(3)+dx]);
end
% Calculate displacement, cost and efficiency of a gait
% Note that, for inertial cost, cost is returned as the integral of torque
% squared, while for drag-based systems, cost is the path length of the
% gait
dcost = zeros(len_x,1);
[~, net_disp_opt, cost] = evaluate_displacement_and_cost1(s,p,[0, T],'interpolated','fixed_step');
[~, dnet_disp_opt1, dcost(1)] = evaluate_displacement_and_cost1(s,dp1,[0, T],'interpolated','fixed_step');
[~, dnet_disp_opt2, dcost(2)] = evaluate_displacement_and_cost1(s,dp2,[0, T],'interpolated','fixed_step');
disp=net_disp_opt(direction); % displacement produced in the chosen direction produced on executing the gait measured in the optimal coordinates 
ddisp = [dnet_disp_opt1(direction);dnet_disp_opt2(direction)];

if len_x == 3
    [~, dnet_disp_opt3, dcost(3)] = evaluate_displacement_and_cost1(s,dp3,[0, T],'interpolated','fixed_step');
    ddisp = [ddisp;dnet_disp_opt3(direction)];
end

%% calcuating ydot so that lagrange equation is zero.
grad = struct();

% Calculate the gradient of 
grad.disp = (ddisp-disp)/dx;
grad.stroke = (dcost-cost)/dx;

% Calculate the lagrange multiplier at the current fourier coefficient.
lambda = pinv(grad.disp)*grad.stroke;

% The ODE is the gradient of the lagrange equation with slight offset
% of the lagrange multiplier.
xdot = grad.stroke-(1-0.02)*lambda*grad.disp;

% The gradient should move toward shirinking the parameter.
xdot = -xdot;
if isnan(xdot)
    xdot = zeros(size(xdot));
end
if ~isreal(xdot)
    1
end
    

end