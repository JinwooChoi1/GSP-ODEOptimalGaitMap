function [A,Aeq]=nonlcon_ellipse_gait(x,s,n,dimension,lb,ub,direction)
%%%%%%%%% 
%
%This function imposes the nonlinear constraint that all the points forming the gait stay within bounds
%
%Inputs:
%
%x: r1 r2 for elliptical axis
%s: System file which contains the connection vector field, CCF's and
%   metric data
%n: Number of points used to parametrize the gaits in a direct
%   transcription method
%dimension: Indicates the number of shape variables of the system
%lb: Lower bound of shape variables for each point which is obtained from the grid inside which an optimal gait is desired
%ub: Upper bound of shape variables for each point which is obtained from the grid inside which an optimal gait is desired
% 
%%%%%%%%%

% % The first step is to obtain a direct transciption parametrization of the gait from the 
% % fourier series parametrization
p= make_ellipse_gait(x);
T = 2*pi;
y1 = [p.phi_def{1}(linspace(0,T,n).') p.phi_def{2}(linspace(0,T,n).')];
y1(end+1,:) = y1(1,:);
y2=y1(:);

%b=length(y2);

% A1 and A2 together impose the constraint that all the points forming the gait stay in bounds
A1=y2+lb;
A2=-y2-ub;

A = [A1;A2];
Aeq = 0;
end
