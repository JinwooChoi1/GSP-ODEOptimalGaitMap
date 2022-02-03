function [value,isterminal,direction] = event_optimal_map_ellipse_gait(x,s,npoints,dimension,direction)
%%%%%%%%%
%
% This function will record each step-optimal gaits and stop ode solver.
%
% Inputs:
%
% x: r1 r2 for elliptical axis
% s: System file which contains the connection vector field, CCF's and
%   metric data
% n: Number of points used to parametrize the gaits in a direct
%   transcription method
% dimension: Indicates the number of shape variables of the system
%   Outputs
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
    
    p = make_ellipse_gait(x);
    [~, temp_disp, ~] = evaluate_displacement_and_cost1(s,p,[0, 2*pi],'interpolated','fixed_step');
    disp = abs(temp_disp(direction));
    
    % The value that we want to be zero
    value = 1;
    % Record the fourier coefficient at step-optimal gaits.
    if abs(disp/bestDisp-3/4) < 1e-2
        stepOptimalGaits{2} = x;
    elseif abs(disp/bestDisp-2/4) < 1e-2
        stepOptimalGaits{3} = x;
    elseif abs(disp/bestDisp-1/4) < 1e-2
        stepOptimalGaits{4} = x;
        value = 0;
    end
    if disp < 1e-6
        value = 0;
    end
    % Halt integration     
    isterminal = 1;
    % The zero can be approached from either direction
    direction = 0;   
end

