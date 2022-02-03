function [value,isterminal,direction] = event_optimal_map(y,s,dimension,direction)
%%%%%%%%%
%
% This function will record each step-optimal gaits and stop ode solver.
%
% Inputs:
%
% y: Matrix containing the Fourier series coefficients
% s: System file which contains the connection vector field, CCF's and
%   metric data
% direction: Direction to optimize travel: 1=x,2=y,3=theta
%
% value(i) : a mathematical expression describing the ith event.
%   An event occurs when value(i) is equal to zero.
% isterminal(i) = 1 if the integration is to terminate
%   when the ith event occurs. Otherwise, it is 0.
% direction(i) = 0 if all zeros are to be located (the default).
%   A value of +1 locates only zeros where the event function is increasing,
%   and -1 locates only zeros where the event function is decreasing.
%   Specify direction = [] to use the default value of 0 for all events.
%%%%%%%%%

global bestDisp stepOptimalGaits;

y = reshape(y,[10 dimension]);
w1=y(end,1);
T = 2*pi/w1;

p = makeGait(y);
[~, temp_disp, cost] = evaluate_displacement_and_cost1(s,p,[0, T],'interpolated','fixed_step');
disp = abs(temp_disp(direction));

% The value that we want to be zero
value = 1;

% Record the fourier coefficient at step-optimal gaits.
if abs(disp/bestDisp-3/4) < 5e-3
    stepOptimalGaits{2,1} = reshape(y,[10 dimension]);
    stepOptimalGaits{2,2} = disp;
    stepOptimalGaits{2,3} = cost;
elseif abs(disp/bestDisp-2/4) < 5e-3
    stepOptimalGaits{3} = reshape(y,[10 dimension]);
    stepOptimalGaits{3,2} = disp;
    stepOptimalGaits{3,3} = cost;
elseif abs(disp/bestDisp-1/4) < 5e-3
    stepOptimalGaits{4} = reshape(y,[10 dimension]);
    stepOptimalGaits{4,2} = disp;
    stepOptimalGaits{4,3} = cost;
    value = 0;
end

% ode solver stops when the gait displacements become zero.
if disp < 1e-3
    value = 0;
end
% Halt integration
isterminal = 1;
% The zero can be approached from either direction
direction = 0;

end

