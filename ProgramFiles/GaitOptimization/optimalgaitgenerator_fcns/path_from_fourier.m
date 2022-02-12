function y = path_from_fourier(f,n,dimension)
% Returns the shape space parametrization of the gait at n points when provided with
% the fourier coefficients f. The gait is returned as a self-closed gait
% (i.e. the first and last rows of the output y are the same point).
% Inputs:
%   f: Fourier coefficients that parametrize the gait.
%   n: Number of points that should compose the gait, less the final
%       self-closed point
%   dimension: Number of shape variables for system

    y = zeros(n+1,dimension);
    % Determine time period based on value of fourier frequency
    w = f(end,1);
    T = 2*pi/w;
    % Create time vector at which to evaluate points of gait
    t = linspace(0,T,n+1);
    % Evaluate the shape-space parametrization of the gait at every time
    % value in t
    for j=1:dimension
        for i=1:1:n+1
            for k = 1:1:size(f,1)-1
                if k == 1
                    y(i,j) = y(i,j) + f(k,j);
                elseif mod(k,2) == 0
                    y(i,j) = y(i,j) + f(k,j)*cos(floor(k/2)*w*t(i));
                else
                    y(i,j) = y(i,j) + f(k,j)*sin(floor(k/2)*w*t(i));
                end
            end
        end
    end
end