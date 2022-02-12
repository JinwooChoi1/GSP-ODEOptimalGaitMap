function [value,isterminal,direction] = event_optimal_map(y,s,dimension,direction,nfparam)
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

global currentDisp;

% The value that we want to be zero
value = 1;

% ode solver stops when the gait displacements become zero.
if currentDisp < 1e-3
    value = 0;
end
% Halt integration
isterminal = 1;
% The zero can be approached from either direction
direction = 0;

end
